#include <Windows.h>
#include <iostream>
#include <GL/glew.h>
#include <GL/GL.h>
#include <GL/freeglut.h>

#define GLFW_INCLUDE_GLU
#define GLFW_DLL
#include <GLFW/glfw3.h>
#include <vector>
#include <cmath>

#define GLM_SWIZZLE
#include <glm/glm.hpp>
#include <glm/gtc/constants.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtx/string_cast.hpp> 

using namespace glm;

// 전역 변수
int Width = 512;
int Height = 512;
std::vector<float> OutputImage;

// 카메라 설정
const vec3 eye(0.0f, 0.0f, 0.0f);
const vec3 u(1.0f, 0.0f, 0.0f);
const vec3 v(0.0f, 1.0f, 0.0f);
const vec3 w(0.0f, 0.0f, 1.0f);
const float l = -0.1f, r = 0.1f, b = -0.1f, t = 0.1f, d = 0.1f;
const int nx = 512, ny = 512;

// Light source (Q1 requirement)
const vec3 lightPos(-4.0f, 4.0f, -3.0f);
const vec3 lightColor(1.0f, 1.0f, 1.0f); // White light
const float lightIntensity = 1.0f;

// 구조체 클래스 정의
struct Ray {
    vec3 origin;
    vec3 direction;
};

struct Material { // Material properties (Q1 requirement)
    vec3 ka;        // Ambient coefficient
    vec3 kd;        // Diffuse coefficient
    vec3 ks;        // Specular coefficient
    float shininess; // Specular power
};

struct Intersection {
    bool hit;
    float distance;
    vec3 point;
    vec3 normal;
    Material material;
};

class Plane {
public:
    vec3 point;
    vec3 normal;
	Material material; // Material properties (Q1 requirement)

    Plane(vec3 p, vec3 n, Material mat) : point(p), normal(n), material(mat) {}

    Intersection intersect(const Ray& ray) const {
        Intersection result = { false, 0, vec3(0), vec3(0), material };
        float denom = dot(normal, ray.direction);
        if (abs(denom) > 1e-6) {
            float t = dot(point - ray.origin, normal) / denom;
            if (t > 0) {
                result.hit = true;
                result.distance = t;
                result.point = ray.origin + t * ray.direction;
				result.normal = normal; // Normal is constant for a plane
            }
        }
        return result;
    }
};

class Sphere {
public:
    vec3 center;
    float radius;
	Material material; // Material properties insted of color

    Sphere(vec3 c, float r, Material mat) : center(c), radius(r), material(mat) {}

    Intersection intersect(const Ray& ray) const {
        Intersection result = { false, 0, vec3(0), vec3(0), material };
        vec3 oc = ray.origin - center;
        float a = dot(ray.direction, ray.direction);
        float b = 2.0f * dot(oc, ray.direction);
        float c = dot(oc, oc) - radius * radius;
        float discriminant = b * b - 4 * a * c;
        if (discriminant > 0) {
            float t = (-b - glm::sqrt(discriminant)) / (2.0f * a);
            if (t > 0) {
                result.hit = true;
                result.distance = t;
                result.point = ray.origin + t * ray.direction;
                result.normal = normalize(result.point - center);
            }
        }
        return result;
    }
};
//장면 클래스
class Scene {
public:
    std::vector<Plane> planes;
    std::vector<Sphere> spheres;

    void addPlane(const Plane& plane) {
        planes.push_back(plane); 
    }

    void addSphere(const Sphere& sphere) {
        spheres.push_back(sphere); 
    }

    Intersection trace(const Ray& ray, float tMin, float tMax) const {
        Intersection closest = { false, tMax, vec3(0), vec3(0), {} };
        for (const auto& plane : planes) {
            Intersection intersection = plane.intersect(ray);
            if (intersection.hit && intersection.distance > tMin && intersection.distance < closest.distance) {
                closest = intersection;
            }
        }
        for (const auto& sphere : spheres) {
            Intersection intersection = sphere.intersect(ray);
            if (intersection.hit && intersection.distance > tMin && intersection.distance < closest.distance) {
                closest = intersection;
            }
        }
        return closest;
    }
};

class Camera {
public:
    vec3 eye;
    vec3 u, v, w;
    float l, r, b, t, d;
    int nx, ny;

    Camera(vec3 e, vec3 u, vec3 v, vec3 w, float l, float r, float b, float t, float d, int nx, int ny)
        : eye(e), u(u), v(v), w(w), l(l), r(r), b(b), t(t), d(d), nx(nx), ny(ny) {
    }

    Ray getRay(int ix, int iy) const {
        float u_s = l + (r - l) * (ix + 0.5f) / nx;
        float v_s = b + (t - b) * (iy + 0.5f) / ny;
        vec3 direction = normalize(u_s * u + v_s * v - d * w);
        return { eye, direction };
    }
};

// Phong Shading Function (Q1 requirement)
vec3 computePhong(const Intersection& intersection, const vec3& eye, const Scene& scene) {
   if (!intersection.hit) return vec3(0.0f); // Black background

   vec3 N = normalize(intersection.normal);
   vec3 L = normalize(lightPos - intersection.point);
   vec3 V = normalize(eye - intersection.point);
   vec3 R = reflect(-L, N);

   // Shadow ray (Q1 requirement)
   Ray shadowRay = { intersection.point + N * 0.001f, L };
   Intersection shadowIntersection = scene.trace(shadowRay, 0.0f, length(lightPos - intersection.point));
   bool inShadow = shadowIntersection.hit;

   // Ambient
   vec3 color = intersection.material.ka;

   if (!inShadow) {
       // Diffuse
       float diff = max(dot(N, L), 0.0f);
       color += intersection.material.kd * lightColor * diff * lightIntensity;

       // Specular
       float spec = pow(max(dot(R, V), 0.0f), intersection.material.shininess);
       color += intersection.material.ks * lightColor * spec * lightIntensity;
   }

   return clamp(color, 0.0f, 1.0f);
}

// Rendering Function (Updated for Q1)
void render(const Scene& scene, const Camera& camera) {
    OutputImage.clear();
    OutputImage.resize(nx * ny * 3);

    for (int j = 0; j < ny; ++j) {
        for (int i = 0; i < nx; ++i) {
            Ray ray = camera.getRay(i, j);
            Intersection inter = scene.trace(ray, 0.0f, std::numeric_limits<float>::max());
            vec3 color = computePhong(inter, camera.eye, scene);

            int idx = (j * nx + i) * 3;
            OutputImage[idx] = color.r;
            OutputImage[idx + 1] = color.g;
            OutputImage[idx + 2] = color.b;
        }
    }
}

// OpenGL Callbacks
void resize_callback(GLFWwindow*, int nw, int nh) {
    Width = nw;
    Height = nh;
    glViewport(0, 0, nw, nh);
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    glOrtho(0.0, static_cast<double>(Width), 0.0, static_cast<double>(Height), 1.0, -1.0);
}

// Main Function
int main(int argc, char* argv[]) {
    // GLFW Initialization
    if (!glfwInit()) {
        std::cerr << "Failed to initialize GLFW" << std::endl;
        return -1;
    }

    GLFWwindow* window = glfwCreateWindow(Width, Height, "Ray Tracer", NULL, NULL);
    if (!window) {
        std::cerr << "Failed to create GLFW window" << std::endl;
        glfwTerminate();
        return -1;
    }

    glfwMakeContextCurrent(window);
    if (glewInit() != GLEW_OK) {
        std::cerr << "Failed to initialize GLEW" << std::endl;
        glfwTerminate();
        return -1;
    }

    glPixelStorei(GL_UNPACK_ALIGNMENT, 1);
    glPixelStorei(GL_PACK_ALIGNMENT, 1);
    glfwSetFramebufferSizeCallback(window, resize_callback);

    // Scene Setup with Q1 Specifications
    Scene scene;

    // Plane P
    Material planeMat{ vec3(0.2f, 0.2f, 0.2f), vec3(1.0f, 1.0f, 1.0f), vec3(0.0f), 0.0f };
    scene.addPlane(Plane(vec3(0.0f, -2.0f, 0.0f), vec3(0.0f, 1.0f, 0.0f), planeMat));

    // Sphere S1
    Material s1Mat{ vec3(0.2f, 0.0f, 0.0f), vec3(1.0f, 0.0f, 0.0f), vec3(0.0f), 0.0f };
    scene.addSphere(Sphere(vec3(-4.0f, 0.0f, -7.0f), 1.0f, s1Mat));

    // Sphere S2
    Material s2Mat{ vec3(0.0f, 0.2f, 0.0f), vec3(0.0f, 0.5f, 0.0f), vec3(0.5f, 0.5f, 0.5f), 32.0f };
    scene.addSphere(Sphere(vec3(0.0f, 0.0f, -7.0f), 2.0f, s2Mat));

    // Sphere S3
    Material s3Mat{ vec3(0.0f, 0.0f, 0.2f), vec3(0.0f, 0.0f, 1.0f), vec3(0.0f), 0.0f };
    scene.addSphere(Sphere(vec3(4.0f, 0.0f, -7.0f), 1.0f, s3Mat));

    // Camera Setup
    Camera camera(eye, u, v, w, l, r, b, t, d, nx, ny);

    // Render Scene
    render(scene, camera);

    // Main Loop
    while (!glfwWindowShouldClose(window)) {
        glClear(GL_COLOR_BUFFER_BIT);
        glDrawPixels(Width, Height, GL_RGB, GL_FLOAT, &OutputImage[0]);
        glfwSwapBuffers(window);
        glfwPollEvents();
        if (glfwGetKey(window, GLFW_KEY_ESCAPE) == GLFW_PRESS || glfwGetKey(window, GLFW_KEY_Q) == GLFW_PRESS) {
            glfwSetWindowShouldClose(window, GL_TRUE);
        }
    }

    glfwDestroyWindow(window);
    glfwTerminate();
    return 0;
}