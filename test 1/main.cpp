#include <string> // For std::string and std::getline
#include <iostream>


int nameAgeSum(int x, int y)
{
    return x + y;
}

int main()
{
    std::cout << "Enter your full name: ";
    std::string name{};
    std::getline(std::cin >> std::ws, name); // read a full line of text into name

    std::cout << "Enter your age: ";
    int age{};
    std::cin >> age;

    std::cout << "Your age + length of name is: " << nameAgeSum(static_cast<int>(name.length()),age) << '\n';

    return 0;
}