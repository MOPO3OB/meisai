Framework modifying pixel data of PNG files. Allows to hide/extract/execute data.

Usage: `meisai [options] image_file <file_to_hide>`

Options:
```
-t, --text TEXT                  Text to hide
-o, --output FILE                Output file
-x, --execute                    Execute the hidden data
-a, --alpha                      Use alpha channel
-c, --compression LEVEL          Compression level [0-9], default: 6
-d, --bit-depth DEPTH            Bit depth [1,2 or 4], default: 1
-i, --info                       Only print info about the file
    --verbose                    Verbose output
-q, --quite                      Hide any output except for result
-h, --help                       Show help message
-v, --version                    Show version
```
