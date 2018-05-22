LICENSE:

This program is free software. It comes without any warranty, to
the extent permitted by applicable law. You can redistribute it
and/or modify it under the terms of the 

GNU GENERAL PUBLIC LICENSE (GPL) Version 3 

that you can find at:

http://159.149.160.51/modtools/downloads/weederlicense.txt

INSTALL

To compile Weeder 2.0 you need a C++ compiler. Usually g++ is present in many Linux distributions or easily installable with the distribution packet manager. 
Mac OS X users can get g++ by installing the free Xcode package from the Apple App Store.  

Check if g++ is available in your system by typing:

g++ -v

Whether the output is not something similar to "Command not found" you should be OK.

To compile Weeder2.0 type:

g++ weeder2.cpp -o weeder2 -O3

An executable binary file named "weeder2" should appear in the same folder.

To work the executable needs the oligo frequency files in the "FreqFiles" subdirectory. 

You can move or copy the executable in any directory you like provided that the FreqFiles folder is moved or copied to the same directory as well.

For a quick help on how to run the program just type:

./weeder2 -h
