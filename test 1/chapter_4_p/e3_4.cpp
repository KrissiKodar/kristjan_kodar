#include <string> // For std::string and std::getline
#include <iostream>

double enterDouble()
{
    std::cout << "Enter a floating point number: ";
    double x{};
    std::cin >> x;
    return x;
}

char enterChar()
{
    std::cout << "Enter one of the following: +, -, * or /: ";
    char x{};
    std::cin >> x;
    return x;
}

void operation(char o, double x, double y)
{
    if (o == '+')
        std::cout << x << " + " << y << " is " << x + y;
    else if (o == '-')
        std::cout << x << " - " << y << " is " << x - y;
    else if (o == '*')
        std::cout << x << " * " << y << " is " << x * y;
    else if (o == '/')
        std::cout << x << " / " << y << " is " << x / y;
}

int main()
{

    double x{enterDouble()};
    double y{enterDouble()};
    char o{enterChar()};

    operation(o, x, y);

    return 0;
}