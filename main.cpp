#include <iostream>
#include <cmath>
#include <stdio.h>
#include <stdlib.h>
#include <vector>
#include <iomanip>

#include <embree3/rtcore.h>
#include <embree3/rtcore_ray.h>

#include "vector3.h"
#include "ray.h"
#include "bxdf.h"
#include "primitive.h"
#include "utils.h"

// Initial scene representation
std::vector<GeometricPrimitive*> scene_geometry;

// Camera settings
struct Camera{
    Vec3 origin;
    Vec3 dir;
    Ray  center_ray;
    Vec3 x_axis;
    Vec3 y_axis ;

    int width;
    int height;
};

// TODO(?): Implement at least stochastic sampling or stratified sampling (maybe something more complex? filtering? quasi-mc?)
// TODO(?): You may also want to implement a proper camera class.
// Generation of the initial sample
inline Vec3 generate_sample(int pixel_x, int pixel_y, int filter_dx, int filter_dy, Camera camera){
    Vec3 d = camera.x_axis * ((float)(pixel_x + get_uniform()-0.5) / camera.width - 0.5) +
             camera.y_axis * ((float)(pixel_y + get_uniform()-0.5) / camera.height - 0.5) + camera.center_ray.dir;
    return normalize(d);
}


// The main part that constructs light paths and actually solves the rendering equation
Vec3 integrate(const RTCScene& embree_scene, const Ray &ray, int max_depth) {

    Ray normal_ray = ray;
    RTCRayHit rtc_ray;

    // The final radiance (the left part of the rendering equation, L_o)
    Vec3 L(0.0);
    // The radiance from the current path segment (L_i in the rendering equation)
    Vec3 F(1.0);

    // Creating intersection context
    RTCIntersectContext context;
    rtcInitIntersectContext(&context);

    for (int depth = 0; depth < max_depth; ++depth) {
        // Setting RTCRay structure from the current Ray
        Ray_to_RTCRayHit(normal_ray, rtc_ray);

        int id = -1;
        // Test it for intersection with the scene
        rtcIntersect1(embree_scene, &context, &rtc_ray);

        // If the ray doesn't hit anything, return the default value (usually Vec3(0.0, 0.0, 0.0)) for the current hit
        // (meaning we stop tracing the path, but it still can carry the energy, depending on your implementation)
        if (rtc_ray.hit.geomID == (int) RTC_INVALID_GEOMETRY_ID) {
            return L + Vec3(0.0);
        }

        // The hit occured and we get the primitive's Id
        id = rtc_ray.hit.geomID;
        // Get the instance from the Id
        const GeometricPrimitive* shape = scene_geometry[id];

        // Since we can only get the distance from the hit, calculate the actual hit point
        Vec3 p = normal_ray(rtc_ray.ray.tfar);

        Vec3 normal;
        Vec3 flipped_normal;
        Vec3 new_dir, bas1, bas2, bas3;
        
        double rnd1, rnd2, chance, cos_theta, sin2_trans, cos_trans, refr, n1 = 1.0, n2 = 1.5, n3, rs, rp;

        // We interpolate normals to get a smooth shading (a.k.a. Gouraud shading)
        // You can set this flag in the mesh's constructor
        if(shape->if_interpolated) {
            float inter_normal[] = {0.0f, 0.0f, 0.0f};
            rtcInterpolate0(rtcGetGeometry(embree_scene, id), rtc_ray.hit.primID, rtc_ray.hit.u, rtc_ray.hit.v,
                            RTC_BUFFER_TYPE_VERTEX_ATTRIBUTE, 0, inter_normal, 3);
            normal = normalize(Vec3f(inter_normal));
            flipped_normal = dot(normal, normal_ray.dir) < 0 ? normal : -normal;
        }
        else {
            normal = normalize(Vec3(rtc_ray.hit.Ng_x, rtc_ray.hit.Ng_y, rtc_ray.hit.Ng_z));
            flipped_normal = dot(normal, normal_ray.dir) < 0 ? normal : -normal;
        }

        // Adding the emmission of the hit primitive (L_e in the rendering equation)
        L += F * shape->emission;
        F *= shape->color;

        // TODO: delete this DONE
        //return shape->color;

        // TODO(?): Implement the Russian roulette here
        if (depth > 4) {
			chance = get_uniform();
			if (chance < 0.2) break;
        }

        // Next path segment
        switch (shape->bxdf) {
            // TODO: Implement specular bounce
            case BxDF_TYPE ::SPECULAR: {
                // Hint: Both PDF/BRDF are 1.0, so we don't calculate them explicitly
                // Hint: shift the hitpoint a little bit along the normal to avoid self-intersection (same with other materials)
                // (something like p += flipped_normal * D_OFFSET_CONSTANT;)
                cos_theta = dot(flipped_normal, -normal_ray.dir);
                new_dir = 2*cos_theta * flipped_normal + normal_ray.dir;
                p += flipped_normal * D_OFFSET_CONSTANT;
                normal_ray = Ray(p, normalize(new_dir));
                break;
            }
            // TODO: Implement refractive bounce
            case BxDF_TYPE::REFRACTIVE: {
                // HINT: Discrete probability of picking the direction and Fresnel reflection are equal in this case
                // and cancel each other, so depending on your implementation you may avoid calculating them explicitly
                // HINT: While implementing the correct refractions, pay attention to normals, you may need to flip them in some cases
                cos_theta = -dot(normal, normal_ray.dir);
                if (cos_theta < 0) {
					cos_theta = -cos_theta;
					normal = -normal;
					n3 = n1; n1 = n2, n2 = n3;
				}
                refr = n1 / n2;
                
                sin2_trans = refr*refr*(1 - cos_theta*cos_theta);
                cos_trans = std::sqrt(1 - sin2_trans);
                if (sin2_trans > 1.) {
					std::cerr << 1 << std::endl;
					break;
				}
				
				rs = (n1*cos_theta - n2*cos_trans) / (n1*cos_theta + n2*cos_trans);
				rs *= rs;
				rp = (n2*cos_theta - n1*cos_trans) / (n2*cos_theta + n1*cos_trans);
				rp *= rp;
				chance = get_uniform();
				if (chance < (rs+rp)/2) {
					// reflect
					new_dir = 2*cos_theta * flipped_normal + normal_ray.dir;
					p += flipped_normal * D_OFFSET_CONSTANT;
				}
				else {
					//refract
					new_dir = refr*normal_ray.dir + (refr*cos_theta - cos_trans)*normal;
					p -= normal * D_OFFSET_CONSTANT;
				}
                normal_ray = Ray(p, normalize(new_dir));
                break;
            }
            // TODO: Implement diffuse bounce
            case BxDF_TYPE::DIFFUSE: {
                // HINT: BRDF * cosine_term and PDF cancel each other
                //F *= F;
                rnd1 = get_uniform();
                rnd2 = get_uniform();
                new_dir = cosine_weighted_hemisphere(rnd1, rnd2);
                new_dir = mybasis(new_dir, flipped_normal);
                p += flipped_normal * D_OFFSET_CONSTANT;
                normal_ray = Ray(p, normalize(new_dir));
                break;
            }
            /*
            // TODO(?): Implement some complex BRDF
            case BxDF_TYPE::DIFFUSE: {
                // HINT: BRDF * cosine_term and PDF cancel each other
                break;
            }
            */
            default: {
                // Unknown BxDF material.
                return L + Vec3(0.0);
            }
        }
    }
    return L;
}

void render(const RTCScene& embree_scene, const Camera& camera, int spp, const int max_path_depth){
    // Image vector
    Vec3* c = new Vec3[camera.width * camera.height];
    
    //Progressive rendering!
	for (int step=1; step <= 20; ++step) {
		spp = 20;
		// Iterate over all the pixels, first height, then width
		for (int y = 0; y < camera.height; y++){

			fprintf(stderr,"\rRendering (%d spp) %5.2f%%",step*spp,100.0*y/(camera.height-1));

			// TODO: You have to parallelize your renderer DONE?
			// Easy option: just use OpenMP / Intel TBB to run the for loop in parallel (don't forget to avoid writing to the same pixel at the same time by different processes)
			// Not so easy option: run ray-intersection workloads in parallel
			#pragma omp parallel for
			for (int x = 0; x < camera.width; ++x) {
				// Getting pixel's index
				int current_idx = (camera.height - y - 1) * camera.width + x;

				Vec3 r = Vec3(0.0);
				for (int s = 0; s < spp; s++) {
					// Generate direction for the initial ray
					Vec3 d = generate_sample(x, y, 0.0, 0.0, camera);
					// Add light path's contribution to the pixel
					r = r + integrate(embree_scene, Ray(camera.center_ray.org, d), max_path_depth) * (1.0 / spp);

				}
				// You might want to clamp the values in the end
				c[current_idx] =  (Vec3(step-1) * c[current_idx] + Vec3(clamp(r.x), clamp(r.y), clamp(r.z))) / Vec3(float(step));
			}
		}
		//~ for (int y = 0; y < camera.height; y++)
			//~ for (int x = 0; x < camera.width; ++x) {
				//~ int current_idx = (camera.height - y - 1) * camera.width + x;
				//~ c[current_idx] = c[current_idx] / Vec3(float(step));
			//~ }
				
		// Save .ppm image.
		// TODO(?): you may want implement your own ldr/hdr image i/o routines
		save_ppm(camera.width, camera.height, c, spp*step);
		
		//~ for (int y = 0; y < camera.height; y++)
			//~ for (int x = 0; x < camera.width; ++x) {
				//~ int current_idx = (camera.height - y - 1) * camera.width + x;
				//~ c[current_idx] = c[current_idx] * Vec3(float(step));
			//~ }
	}
    delete [] c;
}

int main(int argc, char *argv[]){

    std::cout << std::setprecision(16);

    // Image resolution, SPP
    int film_width        = 800;
    int film_height       = 800;
    int samples_per_pixel = 100;

    // Setting the camera
    Camera camera_desc;
    camera_desc.origin        = Vec3(0.0, 0.0, -120.0);
    camera_desc.dir           = Vec3(0.0, 0.0, 1.0);
    camera_desc.center_ray    = Ray(camera_desc.origin, normalize(camera_desc.dir));
    camera_desc.x_axis        = Vec3(1.0, 0.0, 0.0);
    camera_desc.y_axis        = normalize(cross(camera_desc.x_axis, -camera_desc.center_ray.dir));
    camera_desc.width         = film_width;
    camera_desc.height        = film_height;


    // Path to our models
    std::string dragon  = "assets/dragon.obj";
    std::string bunny   = "assets/bunny.obj";

    // Zero transform and transforms for our models
    Transform zero_trans    = Transform(Vec3(1.0, 1.0, 1.0), Vec3(0.0), Vec3(0.0));
    Transform dragon_trans  = Transform(Vec3(70.0, 70.0, 70.0), Vec3(0.0), Vec3(0.0, -30.0, 0.0));
    //Transform bunny_trans   = Transform(Vec3(40.0, 40.0, 40.0), Vec3(0.0), Vec3(10.0, -50.0, 0.0));
    Transform bunny_trans   = Transform(Vec3(40.0, 40.0, 40.0), Vec3(0.0), Vec3(10.0, -40.0, 0.0));

    const double        r = 10000.0;
    const double offset_r = 10050.0;

    // Cornell box
    scene_geometry.push_back(static_cast<GeometricPrimitive*>(new Sphere(r,      zero_trans, Vec3(-offset_r, 0.0, 0.0),         Vec3(00.0), Vec3(0.75,0.25,0.25), DIFFUSE)));//Left
    scene_geometry.push_back(static_cast<GeometricPrimitive*>(new Sphere(r,      zero_trans, Vec3( offset_r, 0.0, 0.0),         Vec3(00.0), Vec3(0.25,0.25,0.75), DIFFUSE)));//Right
    scene_geometry.push_back(static_cast<GeometricPrimitive*>(new Sphere(r,      zero_trans, Vec3(0.0, 0.0, -offset_r - 100.0), Vec3(00.0), Vec3(0.75,0.75,0.75), DIFFUSE)));//Back
    scene_geometry.push_back(static_cast<GeometricPrimitive*>(new Sphere(r,      zero_trans, Vec3(0.0, 0.0,  offset_r),         Vec3(00.0), Vec3(0.75,0.75,0.75),  DIFFUSE)));//Front
    scene_geometry.push_back(static_cast<GeometricPrimitive*>(new Sphere(r,      zero_trans, Vec3(0.0, -offset_r, 0.0),         Vec3(00.0), Vec3(0.75,0.75,0.75), DIFFUSE)));//Bottom
    scene_geometry.push_back(static_cast<GeometricPrimitive*>(new Sphere(r,      zero_trans, Vec3(0.0,  offset_r, 0.0),         Vec3(00.0), Vec3(0.75,0.75,0.75), DIFFUSE)));//Top
    scene_geometry.push_back(static_cast<GeometricPrimitive*>(new Sphere(5000.0, zero_trans, Vec3(0.0, 5049.99, 0.0),           Vec3(12.0), Vec3(0.0),            DIFFUSE)));//Light
    //~ scene_geometry.push_back(static_cast<GeometricPrimitive*>(new Sphere(15.0,    zero_trans, Vec3(0.0, 30, 0.0),                Vec3(12.0), Vec3(0.0),            DIFFUSE)));//Light

    // Other objects
    scene_geometry.push_back(static_cast<GeometricPrimitive*>(new TriangleMesh(dragon, dragon_trans, Vec3(0.0), Vec3(0.0), Vec3(212/255.0, 175/255.0, 55/255.0), DIFFUSE)));

    scene_geometry.push_back(static_cast<GeometricPrimitive*>(new Sphere(10.0,      zero_trans, Vec3(20.0,  10.0, 10.0),         Vec3(00.0), Vec3(1.0, 1.0, 1.0), SPECULAR)));//Top
    scene_geometry.push_back(static_cast<GeometricPrimitive*>(new Sphere(15.0,      zero_trans, Vec3(-30.0,  10.0, 0.0),         Vec3(00.0), Vec3(1.0, 1.0, 1.0), REFRACTIVE)));//Top


    // Creating a new device
    RTCDevice rtc_device = rtcNewDevice("threads=0");

    // Creating a new scene
    RTCScene rtc_scene = rtcNewScene(rtc_device);
    //rtcSetSceneFlags(rtc_scene, RTC_SCENE_FLAG_COMPACT | RTC_SCENE_FLAG_ROBUST);

    // Constructing Embree objects, setting VBOs/IBOs
    for(int i = 0; i < scene_geometry.size(); ++i) {
        scene_geometry[i]->construct_embree_object(rtc_device, rtc_scene);
    }

    // Loading the scene
    rtcCommitScene(rtc_scene);

    // Start rendering
    render(rtc_scene, camera_desc, samples_per_pixel, 10);

    // Releasing the scene and then the device
    rtcReleaseScene(rtc_scene);
    rtcReleaseDevice(rtc_device);

    scene_geometry.clear();
    std::cout << std::endl;
}
