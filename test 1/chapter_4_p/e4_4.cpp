#include <string> // For std::string and std::getline
#include <iostream>

double enterHeight()
{
    std::cout << "Enter then height of the tower in meters: ";
    double x{};
    std::cin >> x;
    return x;
}

double freeFall(double t, double g)
{
    return 0.5*g*t*t;
    
}

void heightOfBall(double x, int t,double g)
{
    double distance{x - freeFall(t,g)}; 
    if (distance > 0)
        std::cout << "At " << t << " seconds, the ball is at height:  " << distance << " meters\n" ;
    else
        std::cout << "At " << t << " seconds, the ball is on the ground.\n";
}

int main()
{
    // simulate ball dropping off tower
    // ask for height of tower in meters
    // ball has no initial velocity
    // output height of ball above ground after 0, 1, 2, 3, 4, and 5 secondsþ
    // ball should not go underneath the ground 
    // formula df=0.5*g*t^2
    const double gravity { 9.8 };
    double initialHeight {enterHeight()};

    heightOfBall(initialHeight, 0, gravity);
    heightOfBall(initialHeight, 1, gravity);
    heightOfBall(initialHeight, 2, gravity);
    heightOfBall(initialHeight, 3, gravity);
    heightOfBall(initialHeight, 4, gravity);
    heightOfBall(initialHeight, 5, gravity);

    return 0;
}