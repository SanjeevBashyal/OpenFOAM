#include "aggregate.H"
#include "point.H"
#include "List.H"
#include "random"

Bashyal::aggregate::aggregate(float r1, float r2, int n)
{
    Foam::List<Foam::point> Nodes(n);
    this.r1 = r1;
    float phi;
    float theta;

    for (int i = 0; i < n; i++)
    {
        phi = 1.0;
        theta = 1.0;
        Foam::Vector<float> pointVector=this.getSurfacePoint(phi,theta)
        Nodes.append()
    }
}

Foam::Vector<float> Bashyal::aggregate::getSurfacePoint(float phi, float theta)
{
    std::random_device rd;                                     // Seed source for the random number engine
    std::mt19937 gen(rd());                                    // Mersenne Twister engine seeded with rd()
    std::uniform_real_distribution<double> realDist(0.0, 1.0); // Real number distribution [0.0, 1.0]
    float x = realDist(gen);
    float y = realDist(gen);
    float z = realDist(gen);
    return Foam::Vector<float>(x, y, z);
}

namespace Bashyal
{
    int hello(int a)
    {
        return a + 1;
    }
}
