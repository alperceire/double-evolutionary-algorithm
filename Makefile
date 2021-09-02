all: clean
	@g++ -o main dea.c -lglut -lGLU -lGL -lm
clean:
	rm main || true
