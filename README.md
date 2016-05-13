

###Compiling:

```bash

    #download maginatics/threadpool to the ./libs directory
    bash scripts/download-install-maginatics-threadpool.sh
    

    mkdir build
    cd build

    #if you are in msys
    GENERATOR="MSYS Makefiles"

    #if you are on a unix-like
    GENERATOR="Unix Makefiles"

    #if you want to compile for msvc, see `cmake --help` for other options
    GENERATOR="Visual Studio 14 2015"

    #For debug mode
    cmake -G"$GENERATOR" .. -DCMAKE_BUILD_TYPE=Debug

    #Or, for release mode
    cmake -G"$GENERATOR" .. -DCMAKE_BUILD_TYPE=Release

    #compile
    cmake --build .


    #later,
    #optionally change the build type to rebuild
    cmake . -DCMAKE_BUILD_TYPE=Debug

    #rebuild
    cmake --build .
```



###Pi hex digits

To generate the pi hex digits:

```bash
    PIHEXDIGITS=100000000
    #PIHEXDIGITS is the number of digits, alter it at will.
    ./gmp-chudnovsky.exe $PIHEXDIGITS 1 16 | cut -c 4- | head --bytes -4

```

Pipe it to a file like so:

```bash

    ./gmp-chudnovsky.exe $PIHEXDIGITS 1 16 | cut -c 4- | head --bytes -4 > pihexdigits

    #check that the number of characters in the file matches PIHEXDIGITS
    wc pihexdigits

```


###Invocation

```bash
    
    #pihexdigits is the previously generated file with pi hex digits
    #input.pizip is some file you want to compress
    #output.pizip is the name of the compressed file you wish to generate
    #chunk size of 4, and 4 threads
    ./pizip -c 4 -t 4 pihexdigits input.text output.pizip

```



