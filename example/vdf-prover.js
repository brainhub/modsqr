const { exec } = require('child_process')

function vdfProver(
    t, // VDF paramter t (to we 2^t squarings)
    x, // input (e.g. sha256 output)
    callback    // results are provided here as the first argument 'result'. 'result' is always 512 hex characters long
) {
    exec('./sqr ' + t + ' ' + x, async (err, out) => {
        if (err) {
            //console.error(err)
            callback("");       // return an empty strings on error. Empty string is never returned on success.
        }
        else {
            callback(out);      // return 2048-bit number, encoded as 512 hexadecimal characters. The 'sqr' always pads to 512 characters. 
        }
    })  
}

t = 20
x = "4adc33bd9fe74303c344be46e5916d65182fb218e248fe80452ab3f025b06c64"

vdfProver(t, x, (result) => {
    console.log( x + "^2^2^" + t + " mod N == " + result );
});

