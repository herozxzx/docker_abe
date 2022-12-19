objs = authorizer.o message_handle.o 

authorizer: $(objs)
	g++ -o ../authorizer $(objs) -L ~/.local/lib -Wl,-rpath ~/.local/lib  -l pbc -lgmp -pthread
authorizer.o: authorizer.cpp ~/.local/include/pbc
	g++ -c authorizer.cpp -I ~/.local/include/pbc -lgmp -pthread
message_handle.o: message_handle.cpp
	g++ -c message_handle.cpp
.PHONY: clean
clean:
	-rm -f $(objs)
