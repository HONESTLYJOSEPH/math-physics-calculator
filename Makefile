# ── Makefile for University Math & Physics Calculator ──────
CC      = gcc
CFLAGS  = -Wall -Wextra -O2 -std=c11
LDFLAGS = -lm
TARGET  = math_physics_calc
SRC     = math_physics_calc.c

.PHONY: all clean run

all: $(TARGET)

$(TARGET): $(SRC)
	$(CC) $(CFLAGS) -o $(TARGET) $(SRC) $(LDFLAGS)
	@echo "Build successful → ./$(TARGET)"

run: all
	./$(TARGET)

clean:
	rm -f $(TARGET)
