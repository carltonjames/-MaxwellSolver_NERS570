#include <iostream>
#include <vector>
#include "Mesh.hpp"
#include "Field.hpp"
#include "Source.hpp"
class Source {
protected:
    float charge;
    int location[3];
    std::vector<float> velocity;
    bool distribution;
    bool dynamic;
    bool free;
    float mass;

public:
    virtual ~Source() = 0;
    virtual float getCharge() const = 0;
    virtual float getCharge(int i, int j, int k) const = 0;
    virtual std::array<float, 3> getLocation() const = 0;
    virtual std::array<float, 3> getVelocity() const = 0;
    virtual bool isDistribution() const = 0;
    virtual bool isDynamic() const = 0;
    virtual bool isFree() const = 0;
    virtual float getMass() const = 0;
    virtual float getMass(int i, int j, int k) const = 0;

    float getCharge() const {
        return charge;
    }

    float getmass() const {
        return mass;
    }

    std::array<float, 3> getLocation() const {
        return location;
    }

    void setVelocity(std::array<float, 3> newVelocity) {
        velocity = newVelocity;
    }

    void setLocation(std::array<float, 3> newLocation) {
        location = newLocation;
    }

    void isFree() {
        return free;
    }

    void isDynamic() {
        return dynamic;
    }

    void isDistribution() {
        return distribution;
    }
};
