#include <iostream>

// Loop between 5 and 1
int main()
{
	int outer{ 1 };
	while (outer <= 5)
	{
		// loop between outer and 1
		int inner{ 5 };
		while (inner >= 1)
        {   
            if (inner > outer)
			    std::cout << "  ";
            else
                std::cout << inner << ' ';
            --inner;
        }

		// print a newline at the end of each row
		std::cout << '\n';
		++outer;
	}

	return 0;
}