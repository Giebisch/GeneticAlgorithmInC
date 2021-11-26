all: build run

build: main.c
	gcc -g main.c -o main -lm
run: main
	./main