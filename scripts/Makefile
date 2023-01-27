dabe = dabe.o message_handle.o
cpabe_option = cpabe_option.o message_handle.o
dabe_option = dabe_option.o message_handle.o
objs = test.o dabe.o cpabe_option.o dabe_option.o message_handle.o

authorizer: $(objs)
	g++ -o ../test test.o -L ~/.local/lib -Wl,-rpath ~/.local/lib  -l pbc -lgmp -pthread
	g++ -o ../dabe $(dabe) -L ~/.local/lib -Wl,-rpath ~/.local/lib  -l pbc -lgmp -pthread
	g++ -o ../cpabe_option $(cpabe_option) -L ~/.local/lib -Wl,-rpath ~/.local/lib  -l pbc -lgmp -pthread
	g++ -o ../dabe_option $(dabe_option) -L ~/.local/lib -Wl,-rpath ~/.local/lib  -l pbc -lgmp -pthread
test.o: test.cpp ~/.local/include/pbc
	g++ -c test.cpp -I ~/.local/include/pbc -lgmp -pthread
dabe.o: dabe.cpp ~/.local/include/pbc
	g++ -c dabe.cpp -I ~/.local/include/pbc -lgmp -pthread
cpabe_option.o: cpabe_option.cpp ~/.local/include/pbc
	g++ -c cpabe_option.cpp -I ~/.local/include/pbc -lgmp -pthread
dabe_option.o: dabe_option.cpp ~/.local/include/pbc
	g++ -c dabe_option.cpp -I ~/.local/include/pbc -lgmp -pthread
message_handle.o: message_handle.cpp
	g++ -c message_handle.cpp
.PHONY: clean
clean:
	-rm -f /pbc/scripts/*.o