#include <iostream>
#include <cmath>



int main()
{
    double y{ 0.0 };
    for (int x{ 0 }; x <= 10; ++x)
        {
            std::cout << "cos af " << y << " er " << std::cos(y) << '\n';
            y += 0.1;
        }  
}