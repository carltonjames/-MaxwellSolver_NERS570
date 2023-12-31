#ifndef SOURCE_H
#define SOURCE_H

#include <array>
#include <string>
#include <Field.hpp>
class Source {
protected:
    float charge = 0.0f;
    std::array<float, 3> location = { 0.0f, 0.0f, 0.0f };
    std::array<float, 3> velocity = { 0.0f, 0.0f, 0.0f };
    bool distribution = false;
    bool dynamic = false;
    bool neutral = false;
    float mass = 0.0f;
    bool boolFree = false;
    std::string name;
    char special = ' ';


public:
    float t0 = 10.0f;
    float spread = 6.0f * 100.0f;

    Source() {}
    virtual ~Source() {}

    virtual std::array<float, 3> getVelocity(std::array<float, 3> loc) const;
    std::array<float, 3> getVelocity() const;
    virtual float getCharge(std::array<float, 3> loc) const;
    void setLocation(std::array<float, 3> loc);
    void setCharge(float q);
    void setVelocity(const std::array<float, 3>& vel);
    float getCharge() const;
    bool isDistribution() const;
    bool isDynamic() const;
    bool isNeutral() const;
    float getMass() const;
    bool isFree() const;
    std::array<float, 3> getLocation() const;
    std::string getName() const;
    virtual float getMass(std::array<float, 3> loc) const;
    virtual Field getField(int T);
};

#endif // SOURCE_H
