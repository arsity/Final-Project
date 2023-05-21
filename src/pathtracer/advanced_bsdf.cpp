#include "bsdf.h"

#include <algorithm>
#include <iostream>
#include <utility>

#include "application/visual_debugger.h"

using std::max;
using std::min;
using std::swap;

namespace CGL {

// Mirror BSDF //

    Vector3D MirrorBSDF::f(const Vector3D wo, const Vector3D wi) {
        return Vector3D();
    }

    Vector3D MirrorBSDF::sample_f(const Vector3D wo, Vector3D *wi, double *pdf) {

        // TODO Assignment 7: Part 1
        // Implement MirrorBSDF

        *pdf = 1.0;
        reflect(wo, wi);
        return reflectance / abs_cos_theta(*wi);
    }

    void MirrorBSDF::render_debugger_node() {
        if (ImGui::TreeNode(this, "Mirror BSDF")) {
            DragDouble3("Reflectance", &reflectance[0], 0.005);
            ImGui::TreePop();
        }
    }

// Microfacet BSDF //

    double MicrofacetBSDF::G(const Vector3D wo, const Vector3D wi) {
        return 1.0 / (1.0 + Lambda(wi) + Lambda(wo));
    }

    double MicrofacetBSDF::D(const Vector3D h) {
        // TODO: proj3-2, part 3
        // Compute Beckmann normal distribution function (NDF) here.
        // You will need the roughness alpha.
        return exp(-1. * pow(tan(acos(h.z)) / alpha, 2)) / (PI * pow(alpha, 2) * pow(h.z, 4));
    }

    Vector3D MicrofacetBSDF::F(const Vector3D wi) {
        // TODO: proj3-2, part 3
        // Compute Fresnel term for reflection on dielectric-conductor interface.
        // You will need both eta and etaK, both of which are Vector3D.
        double cosTheta_i = abs_cos_theta(wi);
        Vector3D A = (eta * eta + k * k);
        Vector3D B = 2. * eta * cosTheta_i;
        double C = cosTheta_i * cosTheta_i;

        Vector3D Rs = (A - B + C) / (A + B + C);
        Vector3D Rp = (A * C - B + 1.) / (A * C + B + 1.);
        Vector3D F = (Rs + Rp) / 2.;

        return F;
    }

    Vector3D MicrofacetBSDF::f(const Vector3D wo, const Vector3D wi) {
        // TODO: proj3-2, part 3
        // Implement microfacet model here.
        if (wo.z > 0 && wi.z > 0) {
            Vector3D h = (wo + wi) / (wo + wi).norm();// half vector
            return F(wi) * G(wo, wi) * D(h) / 4. / wo.z / wi.z;
        }
        else {
            return Vector3D(0.);
        }
    }

    Vector3D MicrofacetBSDF::sample_f(const Vector3D wo, Vector3D *wi, double *pdf) {
        // TODO: proj3-2, part 3
        // *Importance* sample Beckmann normal distribution function (NDF) here.
        // Note: You should fill in the sampled direction *wi and the corresponding *pdf,
        //       and return the sampled BRDF value.

        // Sampling half angle vector in spherial coordinates
        Vector2D r = sampler.get_sample();
        double r1 = r.x;
        double r2 = r.y;

        double theta_h = atan(sqrt(-1. * pow(alpha, 2) * log(1. - r1)));
        double phi_h = 2. * PI * r2;
        Vector3D h(sin(theta_h) * cos(phi_h), sin(theta_h) * sin(phi_h), cos(theta_h));

        // Assign wi based on the sampled half vector
        *wi = (2. * dot(h, wo)) * h - wo;
        wi->normalize();

        // Check sampled wi validity
        if ((*wi).z < 0) {
            *pdf = 0.;
            return Vector3D(0.);

        }

        // Determine pdf of having sampled the final wi
        double p_theta =
                (2. * sin(theta_h) / (pow(alpha, 2) * pow(cos(theta_h), 3))) * exp(-1. * pow(tan(theta_h) / alpha, 2));
        double p_phi = 1. / 2. / PI;
        double p_h = p_theta * p_phi / sin(theta_h);
        *pdf = p_h / 4. / dot(*wi, h);

        return MicrofacetBSDF::f(wo, *wi);
    }

    void MicrofacetBSDF::render_debugger_node() {
        if (ImGui::TreeNode(this, "Micofacet BSDF")) {
            DragDouble3("eta", &eta[0], 0.005);
            DragDouble3("K", &k[0], 0.005);
            DragDouble("alpha", &alpha, 0.005);
            ImGui::TreePop();
        }
    }

// Refraction BSDF //

    Vector3D RefractionBSDF::f(const Vector3D wo, const Vector3D wi) {
        return Vector3D();
    }

    Vector3D RefractionBSDF::sample_f(const Vector3D wo, Vector3D *wi, double *pdf) {
        // TODO Assignment 7: Part 1
        // Implement RefractionBSDF

        *pdf = 1.0;
        auto t = refract(wo, wi, ior);
        if (!t)
            return {};

        auto eta = 1.0 / ior;
        if (wo.z < 0)
            eta = ior;

        return transmittance / abs_cos_theta(*wi) / pow(eta, 2.0);
    }

    void RefractionBSDF::render_debugger_node() {
        if (ImGui::TreeNode(this, "Refraction BSDF")) {
            DragDouble3("Transmittance", &transmittance[0], 0.005);
            DragDouble("ior", &ior, 0.005);
            ImGui::TreePop();
        }
    }

// Glass BSDF //

    Vector3D GlassBSDF::f(const Vector3D wo, const Vector3D wi) {
        return Vector3D();
    }

    Vector3D GlassBSDF::sample_f(const Vector3D wo, Vector3D *wi, double *pdf) {

        // TODO Assignment 7: Part 1
        // Compute Fresnel coefficient and either reflect or refract based on it.

        *pdf = 1.0;
        auto t = refract(wo, wi, ior);
        if (!t) {
            reflect(wo, wi);
            return reflectance / abs_cos_theta(*wi);
        }

        auto eta = 1.0 / ior;
        if (wo.z < 0)
            eta = ior;

        auto R_0 = pow((1.0 - eta) / (1.0 + eta), 2.0);
        auto R = R_0 + (1.0 - R_0) * pow(1.0 - abs_cos_theta(wo), 5.0);

        if (coin_flip(R)) {
            reflect(wo, wi);
            *pdf = R;
            return R * reflectance / abs_cos_theta(*wi);
        }
        else {
            *pdf = 1.0 - R;
            return (1 - R) * transmittance / abs_cos_theta(*wi) / pow(eta, 2.0);
        }
    }

    void GlassBSDF::render_debugger_node() {
        if (ImGui::TreeNode(this, "Refraction BSDF")) {
            DragDouble3("Reflectance", &reflectance[0], 0.005);
            DragDouble3("Transmittance", &transmittance[0], 0.005);
            DragDouble("ior", &ior, 0.005);
            ImGui::TreePop();
        }
    }

    void BSDF::reflect(const Vector3D wo, Vector3D *wi) {

        // TODO Assignment 7: Part 1
        // Implement reflection of wo about normal (0,0,1) and store result in wi.

        // consider wo already in local coordinates
        // return wi in local coordinates
        *wi = Vector3D(-wo.x, -wo.y, wo.z);

    }

    bool BSDF::refract(const Vector3D wo, Vector3D *wi, double ior) {

        // TODO Assignment 7: Part 1
        // Use Snell's Law to refract wo surface and store result ray in wi.
        // Return false if refraction does not occur due to total internal reflection
        // and true otherwise. When dot(wo,n) is positive, then wo corresponds to a
        // ray entering the surface through vacuum.

        *wi = Vector3D();
        auto t = dot(wo, Vector3D{0, 0, 1});

        auto eta = 1.0 / ior;  // default entering case
        if (t < 0)  // existing case
            eta = ior;

        auto internal_term = 1.0 - eta * eta * (1.0 - pow(wo.z, 2.0));
        if (internal_term < 0)
            return false;

        wi->x = -eta * wo.x;
        wi->y = -eta * wo.y;
        wi->z = sqrt(internal_term);

        if (t > 0)
            wi->z = -wi->z;

        return true;
    }

} // namespace CGL
