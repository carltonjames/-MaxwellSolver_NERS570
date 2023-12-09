#include "Source.hpp"

std::array<float, 3> Source::getVelocity(std::array<float, 3> loc) const {
    // Implementation for getVelocity with location parameter
    // For example, return velocity if loc matches the source's location
    return (loc == location) ? velocity : std::array<float, 3>{0.0f, 0.0f, 0.0f};
}

std::array<float, 3> Source::getVelocity() const {
    return velocity;
}

float Source::getCharge(std::array<float, 3> loc) const {
    // Implementation for getCharge with location parameter
    // For example, return charge if loc matches the source's location
    return (loc == location) ? charge : 0.0f;
}

void Source::setLocation(std::array<float, 3> loc) {
    location = loc;
}

void Source::setCharge(float q) {
    charge = q;
}

void Source::setVelocity(const std::array<float, 3>& vel) {
    velocity = vel;
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

float Source::getMass() const {
    return mass;
}

bool Source::isFree() const {
    return boolFree;
}

std::array<float, 3> Source::getLocation() const {
    return location;
}

std::string Source::getName() const {
    return name;
}

float Source::getMass(std::array<float, 3> loc) const {
    // Implementation for getMass with location parameter
    // For example, return mass if loc matches the source's location
    return (loc == location) ? mass : 0.0f;
}
