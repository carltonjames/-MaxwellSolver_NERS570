#include <iostream>
#include <vector>
#include <array>
#include "Mesh.hpp"
#include "Field.hpp"
#include "Source.hpp"
class Source {
protected:
    virtual ~Source() {}
    float charge;
    std::array<int, 3> location = { 0, 0, 0 };
    std::array<float, 3> velocity = { 0.0f, 0.0f, 0.0f };
    bool distribution;
    bool dynamic;
    bool free;
    float mass;

public:
    virtual Source() {};
    void setLocation(int i, int j, int k);
    void setCharge(float q);
    virtual float getCharge() const;
    virtual float getCharge(int i, int j, int k) const = 0;
    virtual std::array<float, 3> velocity(int i, int j, int k) const = 0;
    virtual bool isDistribution() const;
    virtual bool isDynamic() const;
    virtual bool isFree() const;
    virtual float getMass() const;
    virtual float getMass(int i, int j, int k) const = 0;

    void Source::setLocation(int i, int j, int k) {
        location = { i, j, k };
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

    bool Source::isFree() const {
        return free;
    }

    float Source::getMass() const {
        return mass;
    }

    std::array<float, 3> Source::getLocation() {
        return location;
    }
};