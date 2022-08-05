#include <iostream>


int sumTo(int x)
{
    int y{0};
    for (int counter{ 1 }; counter <= x; ++counter)
        {
            y += counter;
        }  
    return y;
}

int main()
{
    std::cout << "Enter an integer: ";
    int x{};
    std::cin >> x;

    std::cout << "The sum from 1 to " << x << " is " << sumTo(x) << '\n';
    std::cout << "The sum from 1 to " << x << " is " << x*(x+1)/2;
}