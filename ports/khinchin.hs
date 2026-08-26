-- Khinchin's constant (OEIS A002210) via the accelerated zeta series.
--
-- Same mathematics as ../khinchin_fast.c:
--
--   ln(2) ln(K0) = -sum_{k=2}^{N-1} ln((k-1)/k) ln((k+1)/k)
--                  + sum_{n>=1} (zeta(2n) - 1 - sum_{k=2}^{N-1} k^(-2n)) / n * h(n),
--
--   h(n) = sum_{j=1}^{2n-1} (-1)^(j+1)/j.
--
-- Runs entirely on Haskell's native Integer — which GHC backs with GMP,
-- making this potentially the fastest of the "native bignum" ports — in
-- fixed point (scale 2^wp), using only GHC boot libraries (base, array).
-- The port carries its own fixed-point kit: pi by Machin's formula,
-- logarithms via the atanh(1/q) series, exp by argument-halving plus
-- Taylor.  The even zeta values come from the positive-term recurrence
-- (n + 1/2) zeta(2n) = sum_j zeta(2j) zeta(2n-2j) — but only below
-- n_direct, the C program's two-region split: above it each term is a
-- literal tail sum with no zeta values, shortening the O(M^2)
-- recurrence (a lazily self-referential array) to O(n_direct^2).
-- Guard bits absorb the drift; nothing here is interval-certified.
--
-- Tested with GHC 9.10.3: byte-identical output to the C program at
-- 1000 digits (0.07 s), digit-exact at 10,000 (54 s serial).
--
-- Build and run:  ghc -O2 khinchin.hs && ./khinchin DIGITS [OUTPUT_FILE]
-- With OUTPUT_FILE the decimal value plus a newline is written there
-- (the same contract as ../khinchin_fast.c); otherwise stdout.

import Data.Array (Array, listArray, (!))
import Data.Bits (shiftL, shiftR)
import Data.List (foldl')
import System.Environment (getArgs)

-- atanh(1/q) = sum_{j>=0} 1/((2j+1) q^(2j+1)); atan(1/q) alternates.
arcRecip :: Bool -> Integer -> Integer -> Integer
arcRecip alternating q one = go 1 p0 p0
  where
    p0 = one `div` q
    q2 = q * q
    go :: Integer -> Integer -> Integer -> Integer
    go j power s
      | power == 0 = s
      | otherwise =
          let power' = power `div` q2
              term = power' `div` (2 * j + 1)
              s' = if alternating && odd j then s - term else s + term
           in go (j + 1) power' s'

-- exp(y) for 0 <= y < 1: halve 16 times, Taylor, square back.
expFixed :: Int -> Integer -> Integer -> Integer
expFixed wp one y = iterate square (taylor 1 one one) !! 16
  where
    arg = y `shiftR` 16
    square x = (x * x) `shiftR` wp
    taylor :: Integer -> Integer -> Integer -> Integer
    taylor j term acc
      | term' == 0 = acc
      | otherwise = taylor (j + 1) term' (acc + term')
      where
        term' = ((term * arg) `shiftR` wp) `div` j

-- zeta(2), zeta(4), ..., zeta(2M) by the positive-term recurrence,
-- memoised in a lazily self-referential array.
zetas :: Int -> Int -> Integer -> Array Int Integer
zetas wp m z1 = zs
  where
    zs = listArray (1, m) (z1 : map zeta [2 .. m])
    zeta n =
      let half = (n - 1) `div` 2
          raw = sum [zs ! j * zs ! (n - j) | j <- [1 .. half]]
          doubled = 2 * (raw `shiftR` wp)
          acc =
            if even n
              then doubled + (((zs ! (n `div` 2)) ^ (2 :: Int)) `shiftR` wp)
              else doubled
       in (2 * acc) `div` fromIntegral (2 * n + 1)

-- K0 scaled by 2^wp, together with wp.
khinchinFixed :: Int -> (Integer, Int)
khinchinFixed digits = (expFixed wp one y, wp)
  where
    prec = ceiling (fromIntegral (digits + 2) * logBase 2 10 :: Double) + 32
    wp = prec + ceiling (sqrt (fromIntegral prec) :: Double) + 64
    one = 1 `shiftL` wp
    bigN = max 3 (floor ((fromIntegral wp :: Double) ** 0.35))
    m = ceiling (fromIntegral wp * log 2 / (2 * log (fromIntegral bigN)) :: Double) + 1
    nDirect =
      min
        (ceiling (fromIntegral wp / (2 * (logBase 2 (fromIntegral bigN) + 3)) :: Double))
        (m + 1)

    -- -ln((k-1)/k) ln((k+1)/k) = 4 atanh(1/(2k-1)) atanh(1/(2k+1)).
    corr =
      sum
        [ (4 * (a * b)) `shiftR` wp
        | k <- [2 .. bigN - 1]
        , let a = arcRecip False (2 * k - 1) one
        , let b = arcRecip False (2 * k + 1) one
        ]

    piFixed = 16 * arcRecip True 5 one - 4 * arcRecip True 239 one
    zs = zetas wp (nDirect - 1) (((piFixed * piFixed) `shiftR` wp) `div` 6)

    -- Main accelerated loop, ascending, incremental power table.
    ks = [2 .. toInteger bigN - 1]
    step (s, h, powers) n =
      let tail' = zs ! n - one - sum powers
          s' = s + ((tail' * h) `div` toInteger n) `shiftR` wp
          h' = h + one `div` (2 * toInteger n + 1) - one `div` (2 * toInteger n)
          powers' = [p `div` (k * k) | (p, k) <- zip powers ks]
       in (s', h', powers')
    (bernSum, hAfter, _) =
      foldl' step (corr, one, [one `div` (k * k) | k <- ks]) [1 .. nDirect - 1]

    -- Direct region (the C program's split): literal tail sums over
    -- k in [N, K], no zeta values; h continues ascending, power table
    -- re-seeded per block of 32 terms.
    directGo s h n
      | n > m = s
      | otherwise =
          let lastN = min (n + 31) m
              bigK =
                max bigN $
                  min
                    (floor (2 ** ((fromIntegral wp + 32) / (2 * fromIntegral n)) :: Double) + 1)
                    (16 * bigN + 64)
              dks = [toInteger bigN .. toInteger bigK]
              pw0 = [one `div` (k ^ (2 * n)) | k <- dks]
              inner (s', h', pw) mm =
                let tail' = sum pw
                    s'' = s' + ((tail' * h') `div` toInteger mm) `shiftR` wp
                    h'' = h' + one `div` (2 * toInteger mm + 1) - one `div` (2 * toInteger mm)
                 in (s'', h'', [p `div` (k * k) | (p, k) <- zip pw dks])
              (s1, h1, _) = foldl' inner (s, h, pw0) [n .. lastN]
           in directGo s1 h1 (lastN + 1)

    mainSum = directGo bernSum hAfter nDirect

    ln2 = 2 * arcRecip False 3 one
    y = (mainSum `shiftL` wp) `div` ln2

-- Decimal string with `digits` digits after the point.
khinchin :: Int -> String
khinchin digits =
  let (value, wp) = khinchinFixed digits
      scaled = (value * 10 ^ digits + (1 `shiftL` (wp - 1))) `shiftR` wp
      text = show scaled
   in take 1 text ++ "." ++ drop 1 text

main :: IO ()
main = do
  args <- getArgs
  let digits = case args of (d : _) -> read d; [] -> 100
      result = khinchin digits
  case args of
    (_ : path : _) -> writeFile path (result ++ "\n")
    _ -> putStrLn result
