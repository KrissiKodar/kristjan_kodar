#include "io.h"
#include <iostream>

int readNumber()
{
    std::cout << "Enter an integer: ";
    int x{};
    std::cin >> x;
    return x;
}

void writeAnswer(int z)
{
    std::cout << "Their sum is: " << z << '\n';
}