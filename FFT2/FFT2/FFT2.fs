namespace FFT2

open System
open System.Numerics

module Dsp =
    let pi = Math.PI (* π *)

    let ( +@ ) (z1 : Complex) (z2 : Complex) = Complex.Add(z1,z2)
    let ( -@ ) (z1 : Complex) (z2 : Complex) = Complex.Subtract(z1,z2)
    let ( *@ ) (z1 : Complex) (z2 : Complex) = Complex.Multiply(z1,z2)

    (* ビット反転 *)
    let bit_reversal_index n p =
        Array.init p (fun k -> (n >>> k) % 2)
        |> Array.toList
        |> List.rev
        |> List.mapi (fun k b -> if b = 1 then (1 <<< k) else 0)
        |> List.fold (+) 0

    let bit_reversal p x = Array.init (Array.length x) (fun k -> x.[bit_reversal_index k p]) 
    
    (* 回転子 *)
    let w p n =
        let t = -2. * pi * (float n) / (float (1 <<< p))
        Complex(cos t, sin t)

    let butterfly (n : int) (m : int) (w : Complex) (x : Complex array) =
        let xn = x.[n]
        let xm = x.[m]
        x.[n] <- xn +@ (w *@ xm)
        x.[m] <- xn -@ (w *@ xm)

        ()

    (* 2^p 点の FFT *)
    let fft (m : int) (x : Complex array) =
        let y = bit_reversal m x
        
        (* butterfly *)
        for p = 0 to m-1 do
            for n = 0 to (1 <<< (m-p-1)) - 1 do
                for r = 0 to (1 <<< p) - 1 do
                    let ha = (1 <<< (p+1))
                    let qu = (1 <<< p)
                    butterfly (ha*n+r) (ha*n+qu+r) (w (p+1) r) y
                done;
            done;
        done;

        y
