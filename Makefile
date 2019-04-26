FLAGS = -static-libgcc -static-libstdc++

OBJECTS = main.o

main.exe: $(OBJECTS)
	g++ -std=c++14 $(FLAGS) -o $@ $(OBJECTS)

main.o: main.cpp
	g++ -std=c++14 $(FLAGS) -c -o $@ $<

clean:
	-rm $(OBJECTS) main.exe
