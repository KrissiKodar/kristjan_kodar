#include <string> // For std::string and std::getline
#include <iostream>


bool isEven(int x)
{
    return x%2 == 0;
}

int main()
{

    std::cout << "Enter an integer: ";
    int x{};
    std::cin >> x;

    if (isEven(x))
        std::cout << x << " is even\n";
    else
        std::cout << x << " is odd\n";

    return 0;
}