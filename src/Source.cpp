#include <iostream>
#include <vector>
#include <array>
#include <string>

class Source {
protected:
    float charge = 0.0f;
    std::array<float, 3> location = { 0.0f, 0.0f, 0.0f };
    std::array<float, 3> velocity = { 0.0f, 0.0f, 0.0f };
    bool distribution = false;
    bool dynamic = false;
    bool isFree = false;
    bool neutral = false;
    float mass = 0.0f;
    

public:
    Source() {}
    virtual ~Source() {}

    virtual std::array<float, 3> getVelocity() const = 0;
    virtual float getCharge(std::array<float, 3>) const = 0;
    virtual std::string getName() const = 0;
    virtual float getMass(std::array<float, 3>) const = 0;
    virtual string name;

    void setLocation(std::array<float, 3> loc);
    void setCharge(float q);
    float getCharge() const;
    bool isDistribution() const;
    bool isDynamic() const;
    bool isNeutral() const;
    bool isFree() const;
    float getMass() const;
    std::array<float, 3> getLocation() const;
};

void Source::setLocation(std::array<float, 3> loc) {
    location = loc;
}

void Source::setCharge(float q) {
    charge = q;
}

float Source::getCharge() const {
    return charge;
}

bool Source::isDistribution() const {
    return distribution;
}

bool Source::isDynamic() const {
    return dynamic;
}

bool Source::isNeutral() const {
    return neutral;
}

bool Source::isFree() const {
    return isFree;
}

float Source::getMass() const {
    return mass;
}

std::array<float, 3> Source::getLocation() const {
    return location;
}
