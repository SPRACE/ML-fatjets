g++ generateBkg.cc ../lib/libpythia8.a -o generateBkg -I./ -I../include -O2 -ansi -pedantic -W -Wall -Wshadow -fPIC -Wl,-rpath,../lib -ldl -lz -L./ -Wl,-rpath,./ -lfastjet -L/Users/trtomei/root/lib -lCore -lRIO -lHist -pthread -std=c++11 -m64 -I/Users/trtomei/root/include -DGZIPSUPPORT

g++ generateSig.cc ../lib/libpythia8.a -o generateSig -I./ -I../include -O2 -ansi -pedantic -W -Wall -Wshadow -fPIC -Wl,-rpath,../lib -ldl -lz -L./ -Wl,-rpath,./ -lfastjet -L/Users/trtomei/root/lib -lCore -lRIO -lHist -pthread -std=c++11 -m64 -I/Users/trtomei/root/include -DGZIPSUPPORT
