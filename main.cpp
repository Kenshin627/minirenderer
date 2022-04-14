
#include "tgaimage.h"
#include "geometry.h"
#include "model.h"

void line(Vec2i t0, Vec2i t1, TGAImage& image, TGAColor color);

void scanLine_triangle(Vec2i t0, Vec2i t1, Vec2i t2, TGAImage& image, TGAColor color);

void triangle(Vec3f points[3], TGAImage &image, float n1, float n2, float n3, float *zbuffer);

Vec3f barycentric(Vec3f pts[3], Vec3f p);

void rasterize(Vec2i p0, Vec2i p1, TGAImage& image, TGAColor color, int ybuffer[]);

const TGAColor white = TGAColor(255, 255, 255, 255);

const TGAColor red = TGAColor(255, 0, 0, 255);

Model* model = NULL;
const int width = 1000;
const int height = 1000;

int main(int argc, char** argv) {
#pragma region bresenham line
    /*if (2 == argc) {
            model = new Model(argv[1]);
        }
        else {
            model = new Model("obj/african_head.obj");
        }

        TGAImage image(width, height, TGAImage::RGB);
        for (int i = 0; i < model->nfaces(); i++) {
            std::vector<int> face = model->face(i);
            for (int j = 0; j < 3; j++) {
                Vec3f v0 = model->vert(face[j]);
                Vec3f v1 = model->vert(face[(j + 1) % 3]);
                int x0 = (v0.x + 1.) * width / 2.;
                int y0 = (v0.y + 1.) * height / 2.;
                int x1 = (v1.x + 1.) * width / 2.;
                int y1 = (v1.y + 1.) * height / 2.;
                line(x0, y0, x1, y1, image, white);
            }
        }

        image.flip_vertically();
        image.write_tga_file("output.tga");
        delete model;
        return 0;*/
#pragma endregion

#pragma region scanLine_triangle
    //TGAImage image(width, height, TGAImage::RGB);
    //TGAColor red = TGAColor(255, 0, 0, 255);
    //TGAColor white = TGAColor(255, 255, 255, 255);
    //TGAColor green = TGAColor(0, 255, 0, 255);
    //Vec2i t0[3] = { Vec2i(10, 70),   Vec2i(50, 160),  Vec2i(70, 80) };
    //Vec2i t1[3] = { Vec2i(180, 50),  Vec2i(150, 1),   Vec2i(70, 180) };
    //Vec2i t2[3] = { Vec2i(180, 150), Vec2i(120, 160), Vec2i(130, 180) };
    //scanLine_triangle(t0[0], t0[1], t0[2], image, red);
    //scanLine_triangle(t1[0], t1[1], t1[2], image, green);
    //scanLine_triangle(t2[0], t2[1], t2[2], image, white);
    //image.flip_vertically();
    //image.write_tga_file("output.tga");
    //return 0;
#pragma endregion

#pragma region triangle
    //TGAImage image(width, height, TGAImage::RGB);
    //TGAColor red = TGAColor(255, 0, 0, 255);
    //Vec2i t0[3] = { Vec2i(10, 70),   Vec2i(50, 160),  Vec2i(70, 80) };
    //triangle(t0, image, red);
    //image.flip_vertically();
    //image.write_tga_file("output.tga");
    //return 0;
#pragma endregion

#pragma region back-face-culling
    //Vec3f light_dir(0, 0, -1);
    //TGAImage image(width, height, TGAImage::RGB);
    //model = new Model("obj/african_head.obj");
    //for (int i = 0; i < model->nfaces(); i++)
    //{
    //    std::vector<int> face = model->face(i);
    //    Vec2i screen_coords[3];
    //    Vec3f world_coords[3];
    //    for (int j = 0; j < 3; j++)
    //    {
    //        Vec3f v = model->vert(face[j]);
    //        screen_coords[j] = Vec2i((v.x + 1.) * width / 2., (v.y + 1.) * height / 2.);
    //        world_coords[j] = v;
    //    }

    //    Vec3f n = (world_coords[2] - world_coords[0]) ^ (world_coords[1] - world_coords[0]);
    //    n.normalize();
    //    float intensity = n * light_dir;
    //    if (intensity > 0)
    //    {
    //        triangle(screen_coords, image, TGAColor(intensity * 255, intensity * 255, intensity * 255, 255));
    //    }

    //}
    //image.flip_vertically();
    //image.write_tga_file("output.tga");
    //delete model;
    //return 0;
#pragma endregion

#pragma region y-buffer
    //TGAImage scene(width, 16, TGAImage::RGB);
    //TGAColor red = TGAColor(255, 0, 0, 255);
    //TGAColor blue = TGAColor(0, 0, 255, 255);
    //TGAColor green = TGAColor(0, 255, 0, 255);
    //int ybuffer[width];
    //for (int i = 0; i < width; i++)
    //{
    //    ybuffer[i] = std::numeric_limits<int>::min();
    //}
    //rasterize(Vec2i(20, 34), Vec2i(744, 400), scene, red, ybuffer);
    //rasterize(Vec2i(120, 434), Vec2i(444, 400), scene, green, ybuffer);
    //rasterize(Vec2i(330, 463), Vec2i(594, 200), scene, blue, ybuffer);
    //// screen line
    ////line(Vec2i(10, 10), Vec2i(790, 10), scene, white);

    //scene.flip_vertically(); // i want to have the origin at the left bottom corner of the image
    //scene.write_tga_file("scene.tga");
#pragma endregion

#pragma region z-buffer
    float* zbuffer = new float[width * height];
    for (int i = width * height; i--; zbuffer[i] = -std::numeric_limits<float>::max());
    Vec3f light_dir(0, 0, -1);
    TGAImage image(width, height, TGAImage::RGB);
    model = new Model("obj/african_head.obj");
    for (int i = 0; i < model->nfaces(); i++)
    {
        std::vector<int> face = model->face(i);
        Vec3f screen_coords[3];
        Vec3f world_coords[3];
        float intensity[3];
        for (int j = 0; j < 3; j++)
        {
            Vec3f v = model->vert(face[j]);
            screen_coords[j] = Vec3f((v.x + 1.) * width / 2., (v.y + 1.) * height / 2., v.z);
            world_coords[j] = v;
            intensity[j] = model->norm(i, j) * light_dir;
        }

        Vec3f n = (world_coords[2] - world_coords[0]) ^ (world_coords[1] - world_coords[0]);
        n.normalize();
        //float intensity = n * light_dir;
        
        if (intensity > 0)
        {
            triangle(screen_coords, image, intensity[0], intensity[1], intensity[2], zbuffer);
        }

    }
    image.flip_vertically();
    image.write_tga_file("output.tga");
    delete model;
    return 0;
#pragma endregion


}

void line(Vec2i t0, Vec2i t1, TGAImage& image, TGAColor color) {
    int x0 = t0.x;
    int y0 = t0.y;
    int x1 = t1.x;
    int y1 = t1.y;
    bool steep = false;
    if (std::abs(x1 - x0) < std::abs(y1 - y0))
    {
        std::swap(x0, y0);
        std::swap(x1, y1);
        steep = true;
    }

    if ( x0 > x1)
    {
        std::swap(x0, x1);
        std::swap(y0, y1);
    }
    int dy = y1 - y0;
    int dx = x1 - x0;
    int derror = std::abs(dy) * 2;
    int error = 0;
    int y = y0;
    for (int x = x0; x <= x1; x++)
    {
        if (steep)
        {
            image.set(y, x, color);
        }
        else {
            image.set(x, y, color);
        }
        error += derror;
        if (error > dx)
        {
            y += (y1 > y0 ? 1 : -1);
            error -= dx * 2;
        }
    }
}

void scanLine_triangle(Vec2i t0, Vec2i t1, Vec2i t2, TGAImage& image, TGAColor color) {
    if (t0.y > t1.y)
    {
        std::swap(t0, t1);
    };
    if (t0.y > t2.y)
    {
        std::swap(t0, t2);
    };
    if (t1.y > t2.y)
    {
        std::swap(t1, t2);
    }
    int total_height = t2.y - t0.y;
    for (int y = t0.y; y <= t1.y; y++)
    {
        int segment_height = t1.y - t0.y + 1;
        float alpha = (float)(y - t0.y) / total_height;
        float beta = (float)(y - t0.y) / segment_height;
        Vec2i A = t0 + (t2 - t0) * alpha;
        Vec2i B = t0 + (t1 - t0) * beta;
        if (A.x > B.x)
        {
            std::swap(A, B);
        }
        for (int j = A.x; j < B.x; j++)
        {
            image.set(j, y, color);
        }
    }
    for (int y = t1.y; y <= t2.y; y++)
    {
        int segment_height = t2.y - t1.y;
        float alpha = (float)(y - t0.y) / total_height;
        float beta = (float)(y - t1.y) / segment_height;
        Vec2i A = t0 + (t2 - t0) * alpha;
        Vec2i B = t1 + (t2 - t1) * beta;
        if (A.x > B.x)
        {
            std::swap(A, B);
        }
        for (int j = A.x; j < B.x; j++)
        {
            image.set(j, y, color);
        }
    }

}

void triangle(Vec3f points[3], TGAImage &image, float n1, float n2, float n3, float *zbuffer) {
    Vec2f bboxmin(image.get_width() - 1, image.get_height() - 1);
    Vec2f bboxmax(0, 0);
    Vec2f clamp(image.get_width() - 1, image.get_height() - 1);
    for (int i = 0; i < 3; i++)
    {
        Vec3f p = points[i];
        bboxmin.x = std::max(0.0f, std::min(bboxmin.x, p.x));
        bboxmin.y = std::max(0.0f, std::min(bboxmin.y, p.y));

        bboxmax.x = std::min(clamp.x, std::max(bboxmax.x, p.x));
        bboxmax.y = std::min(clamp.y, std::max(bboxmax.y, p.y));
    }

    for (int i = bboxmin.x; i <= bboxmax.x; i++)
    {
        for (int j = bboxmin.y; j <= bboxmax.y; j++)
        {
            Vec3f p;
            p.x = i;
            p.y = j;
            Vec3f _bartcentric = barycentric(points, p);
            if (_bartcentric.x < 0 || _bartcentric.y <0 || _bartcentric.z < 0)
            {
                continue;
            }
            p.z = 0;
            p.z += points[0].z * _bartcentric.x;
            p.z += points[1].z * _bartcentric.y;
            p.z += points[2].z * _bartcentric.z;
            float intensity = _bartcentric.x * n1 + _bartcentric.y * n2 + _bartcentric.z * n3;
            if (zbuffer[int(p.x+p.y*width)] < p.z)
            {
                zbuffer[int(p.x + p.y * width)] = p.z;
                image.set(p.x, p.y, TGAColor(intensity*255, intensity*255, intensity*255,255));
            }
        }
    }
}

Vec3f barycentric(Vec3f pts[3], Vec3f p) {
    Vec3f u = Vec3f(pts[2].x - pts[0].x, pts[1].x - pts[0].x, pts[0].x - p.x) ^ Vec3f(pts[2].y - pts[0].y, pts[1].y - pts[0].y, pts[0].y - p.y);
    if (std::abs(u.z) < 1) return Vec3f(-1, 1, 1);
    return Vec3f(1.f - (u.x + u.y) / u.z, u.y / u.z, u.x / u.z);
}

void rasterize(Vec2i p0, Vec2i p1, TGAImage& image, TGAColor color, int ybuffer[]) {
    if (p0.x > p1.x)
    {
        std::swap(p0, p1);
    }
    for (int x = p0.x; x <= p1.x; x++)
    {
        float t = (x - p0.x) / (float)(p1.x - p0.x);
        int y = p0.y * (1. - t) + p1.y * t;
        if (ybuffer[x] < y)
        {
            ybuffer[x] = y;
            image.set(x, 0, color);
        }
    }
}