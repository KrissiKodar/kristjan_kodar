#include <iostream>

int calculate(int x, int y, char op)
{
    int z {};

    switch (op)
    {
    case '+':
        z = x + y;
        break;
    case '-':
        z = x - y;
        break;
    case '*':
        z = x * y;
        break;
    case '/':
        z = x / y;
        break;
    case '%':
        z = x % y;
        break;
    default:
        std::cerr << "ERROR INVALID OPERATOR!\n";
        return 0;
    }
    std::cout << x << " " << op << " " << y <<  " = " << z;
    return 0;
}

int main()
{
    std::cout << "Enter an integer: ";
    int x {};
    std::cin >> x;

    std::cout << "Enter another integer: ";
    int y {};
    std::cin >> y;
    
    std::cout << "Enter an operator: +, -, *, / or % (modulus):  ";
    char op {};
    std::cin >> op;

	calculate(x, y, op);

}