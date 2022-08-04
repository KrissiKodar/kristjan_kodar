#include <iostream>

int main()
{
    char c{ 'a' };

    while (c < 'z') // infinite loop
    {
        std::cout << "Letter: " << c << "\t ASCII code: " << static_cast<int>(c) << '\n';
        ++c;
    }

    return 0;
}