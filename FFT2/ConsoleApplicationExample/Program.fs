open System.Numerics
open FFT2.Dsp

[<EntryPoint>]
let main argv =
    let pi = 4. * atan 1.
    let x = Array.init 65536 (fun t -> Complex(cos (2. * pi * (float t) / 65536.0), 0.))
    let y = fft 16 x
    printfn "%f" y.[1].Magnitude

    0 // return an integer exit code
