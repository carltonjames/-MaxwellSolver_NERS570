#ifndef SOURCE_H
#define SOURCE_H

#include <vector>
#include <array>

class Source {
public:
    virtual ~Source() = 0;
    float getCharge() const;
    virtual float getCharge(int i, int j, int k) const = 0;
    std::array<float, 3> getLocation() const;
    virtual std::array<float, 3> getVelocity(int i, int j, int k) const = 0;
    void setLocation(std::array<float, 3> newLocation);
    void setVelocity(std::array<float, 3> newVelocity);
    bool isDistribution() const;
    bool isDynamic() const;
    bool isFree() const;
    float getMass() const;
    virtual float getMass(int i, int j, int k) const = 0;

protected:
    float charge;
    std::array<int, 3> location = { 0, 0, 0 };
    std::array<float, 3> velocity = { 0, 0, 0 };
    bool distribution;
    bool dynamic;
    bool free;
    float mass;
};

#endif // SOURCE_H
