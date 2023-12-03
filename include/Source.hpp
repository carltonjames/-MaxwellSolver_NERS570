#ifndef SOURCE_H
#define SOURCE_H

#include <vector>
#include <array>

class Source {
public:
    virtual ~Source() = 0;
    virtual float getCharge() const = 0;
    virtual float getCharge(int i, int j, int k) const = 0;
    virtual std::array<int, 3> getLocation() const = 0;
    virtual std::array<float, 3> getVelocity(int i, int j, int k) const = 0;
    virtual bool isDistribution() const = 0;
    virtual bool isDynamic() const = 0;
    virtual bool isFree() const = 0;
    virtual float getMass() const = 0;
    virtual float getMass(int i, int j, int k) const = 0;


protected:
    float charge;
    std::array<int, 3> location = { 0, 0, 0 };
    std::array<float, 3> velocity;
    bool distribution;
    bool dynamic;
    bool free;
    float mass;
};

#endif // SOURCE_H
