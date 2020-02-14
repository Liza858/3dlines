#include <iostream>
#include <vector>
#include <string>
#include <fstream>
#include <cstdlib>
#include <math.h>
#include <map>
#include <utility>
#include <opencv/cv.h>
#include <opencv2/highgui/highgui.hpp>
#include <opencv2/imgproc/imgproc.hpp>
#include <opencv2/core/core.hpp>

using std::vector;
using std::string;
using std::ifstream;
using std::ofstream;
using std::cin;
using std::cout;
using std::endl;
using std::map;
using std::pair;

#define MAX_LINE_SIZE 20
#define COUNT_OF_ITERATION 50
#define DIRECTION_ERROR 0.1

class camera {


    private:

        vector<vector<double>> R;
        vector<double> T;
        double f;
        int w;
        int h;
        cv::Mat image;
        vector<vector<vector<double>>> img_lines;
        void detect_image(int number) {
            if (number != 1) {
                cv::Mat image_black(image.rows, image.cols ,CV_32S, cv::Scalar(0,0,0));
		        cv::Ptr<cv::LineSegmentDetector> lsd = cv::createLineSegmentDetector(cv::LSD_REFINE_ADV);
		        vector<cv::Vec4f> lines;
                lsd->detect(image, lines);
                cout << lines.size() << endl;
                for (int i = 0; i < lines.size(); i++){
                    cv::Vec4i line = lines[i];
                    cv::line(image_black,
                    cv::Point(line[0], line[1]),
                    cv::Point(line[2], line[3]),  cv::Scalar(i+1, 0, 0), 1);
                }    
                image = image_black;
             }

            for (int  i = 0; i < h; i++) {
                for (int j = 0; j < w; j++) {
                    int color = static_cast<int>(image.at<int32_t>(i, j));
                    if (color != 0) {
                        img_lines[color].push_back({(j - w/2)*1.0, (h/2 - i)*1.0});
                    }
                }
            }
        }

        vector<double> dot(vector<vector<double>> matrix, vector<double> vec) {
            if (vec.size() != matrix[0].size()) {
                throw "dot : error vectors size!";
            }
            vector<double> out;
            for (int i = 0; i < 3; i++) 
	        {
		        double sum = 0;
		        for (int j = 0; j < 3; j++)
		        {
			        sum += matrix[i][j] * vec[j];
                }
                out.push_back(sum);
	        }
            return out;
        }

        vector<double> add(vector<double> vec1, vector<double> vec2) {
            if (vec1.size() != vec2.size()) {
                throw "add : error vectors size!";
            }
            for (int i = 0; i < vec1.size(); i++) 
	        {
		        vec1[i] += vec2[i];
	        }
            return vec1;
        }

        vector<double> sub(vector<double> vec1, vector<double> vec2) {
            if (vec1.size() != vec2.size()) {
                throw "sub : error vectors size!";
            }
            for (int i = 0; i < vec1.size(); i++) 
	        {
		        vec1[i] -= vec2[i];
	        }
            return vec1;
        }

        vector<vector<double>> inverse(vector<vector<double>> input_matrix) {
            vector<vector<double>> out{{0,0,0}, {0,0,0}, {0,0,0}};
            cv::Mat matrix(3, 3, CV_64F);
            for (int i = 0; i < 3 ; i++) 
            {
                for (int j = 0; j < 3; j++) 
                {
                    matrix.at<double>(i, j) = input_matrix[i][j];
                }
            }
            matrix = matrix.inv();
            for (int i = 0; i < 3; i++) 
	        {
		        for (int j = 0; j < 3; j++)
		        {
			        out[i][j] = matrix.at<double>(i, j);
                }
	        }
            return out;
        }
        
      
         bool check_error(vector<double> dir1, vector<double> dir2, double& error) {
            vector<double> z = sub(dir1, dir2);
            double abss = sqrt(z[0]*z[0] + z[1]*z[1]);
            // cout << "abs1 " << abss << endl;
            if (abss < DIRECTION_ERROR) {
                error = abss;
                return true;
            } else {
                z = add(dir1, dir2);
            }
            abss = sqrt(z[0]*z[0] + z[1]*z[1]);
            // cout << "abs2 " << abss << endl;
            if (abss < DIRECTION_ERROR) {
                error = abss;
                return true;
            }
            error = abss;
            return false;
         }

    public:

        int num;

        camera() 
        {}

        camera(int number, string cameras_filename, string image_file) 
            : num(number)
        {
            image = cv::imread(image_file, CV_LOAD_IMAGE_GRAYSCALE);
            if (number == 1) {
                image.convertTo(image, CV_32S);
            }
            w = image.cols;
            h = image.rows;
            img_lines = vector<vector<vector<double>>>(w*h);
            detect_image(number);
            ifstream in(cameras_filename, std::ifstream::in);
            string line;
            getline(in, line);
            int count_cameras = 0;
            int other = 0;
            in >> count_cameras >> other;
            if (number > count_cameras || number < 1) {
                cout << "error camera number!" << endl;
                exit(1);
            }
            for (int i = 5; i < 5 * number + 1; i++) {
                 getline(in, line);
            }
            double x = 0;
            double y = 0;
            double z = 0;
            in >> x >> y >> z;
            f = x;
            R.clear();
            for (int i = 0; i < 3; i++) {
                 in >> x >> y >> z;
                 vector<double> v = {x, y, z};
                 R.push_back(v);
            } 
            in >> x >> y >> z;
            T = vector<double>{x, y, z};
        }
        ~camera() {}
        vector<vector<double>> get_R() { return R; }
        vector<double> get_T() { return T; }
        double get_f() { return f; }
        int get_w() { return w; }
        int get_h() { return h; }
        vector<vector<vector<double>>>& get_lines() {return img_lines;}
        cv::Mat get_image() { return image; } 
        void print_R() {
            for (int i = 0; i < 3; i++) {
                for (int j = 0; j < 3; j++) {
                    cout << R[i][j] << " ";
                }
                cout << endl;
            }
        }

        void print_T() {
            for (int i = 0; i < 3; i++) {
                cout << T[i] << " ";
            }
            cout << endl;
        }

        vector<double> change_global_coords_to_local(vector<double> global_xyz) {
            vector<double> local_xyz;
            local_xyz = dot(R, global_xyz);
            local_xyz = add(local_xyz, T);
            return local_xyz;
        }
         
        vector<double> change_global_coords_to_pixel(vector<double> global_xyz) {
            vector<double> local_xyz = change_global_coords_to_local(global_xyz);
            if (local_xyz[2] == 0) {
                throw "error!";
            }
            vector<double> pixel_xy = { -local_xyz[0] * f / local_xyz[2], -local_xyz[1] * f / local_xyz[2] };
            return pixel_xy;
        }

        vector<double> change_pixel_to_global_coords(vector<double> pixel_xy, double local_z) {
            vector<double> local_xyz = { -pixel_xy[0] * local_z / f,  -pixel_xy[1] * local_z / f, local_z };
	        vector<vector<double>>R_inv = inverse(R);
	        vector<double> global_xyz = dot(R_inv , sub(local_xyz, T));
	        return global_xyz;
        }

        vector<double> find_plane(vector<double> pixel1, vector<double> pixel2) {
            vector<double> xyz1 = change_pixel_to_global_coords(pixel1, 8);
	        vector<double> xyz2 = change_pixel_to_global_coords(pixel1, 10);
	        vector<double> xyz3 = change_pixel_to_global_coords(pixel2, 12);
            vector<double> x1 = sub(xyz2, xyz1);
            vector<double> x2 = sub(xyz3, xyz1);       
	        vector<double> n = { x1[1]*x2[2] - x1[2]*x2[1], x1[2]*x2[0] - x1[0]*x2[2], x1[0]*x2[1] - x1[1]*x2[0] };
	        double A = n[0];
	        double B = n[1];
	        double C = n[2];
	        double D = - A * xyz2[0] - B * xyz2[1] - C * xyz2[2];
	        return vector<double>{A, B, C, D};
        }

        vector<double> find_line(vector<double> pixel_xy) {
            double z1 = 8.0;
	        double z2 = 5.0;
	        vector<double> xyz1 = change_pixel_to_global_coords(pixel_xy, z1);
	        vector<double> xyz2 = change_pixel_to_global_coords(pixel_xy, z2);
	        return vector<double>{xyz2[0] - xyz1[0], xyz2[1] - xyz1[1], xyz2[2] - xyz1[2], xyz1[0], xyz1[1], xyz1[2]};
        }

        vector<double> intersection(vector<double> plane, vector<double> line) {
            double A = plane[0];
            double B = plane[1];
            double C = plane[2];
            double D = plane[3];
            double a = line[0];
            double b = line[1];
            double c = line[2];
            double x0 = line[3];
            double y0 = line[4];
            double z0 = line[5];
            double func1 = x0 / a * (B * b + C * c) - B * y0 - C * z0 - D;
	        double func2 = A + B * b / a + C * c / a;
	        double x = func1/func2;
	        double y = y0 + b / a * (x - x0);
	        double z = z0 + c / a * (x - x0);
	        return vector<double>{x, y, z};
        }

        vector<double> find_line_by_two_points(vector<double> pixel1, vector<double> pixel2) {
            return vector<double>{pixel2[0] - pixel1[0], pixel2[1] - pixel1[1], pixel1[0], pixel1[1]};

        }

        void print(vector<double> vec) {
            for (int i = 0; i < vec.size(); i++) {
                cout << vec[i] << " ";
            } 
            cout << endl;
        }

        void find_epipolar_line(vector<double> global_xyz1, vector<double> global_xyz2, vector<vector<double>>& out) {
            vector<double> pix1 = change_global_coords_to_pixel(global_xyz1);
            vector<double> pix2 = change_global_coords_to_pixel(global_xyz2);
            vector<double> line = find_line_by_two_points(pix1, pix2);
            int min_pixel_x = - w / 2;
            int max_pixel_x = w / 2;
            int min_pixel_y = - h / 2;
            int max_pixel_y = h / 2;
            for (int i = min_pixel_x; i <= max_pixel_x; i++) {
                double j = 1.0 * line[1] / line[0] * (i - line[2]) + line[3];
                if (j >= min_pixel_y && j <= max_pixel_y) {
                    out.push_back(vector<double>{i*1.0, j});
                }
            }
        }
        
        vector<double> find_pixel(int x, int y, int col) {
            vector<double> pixel{-w - 100.0, 0};
            for (int i = -4; i <= 4; i++ ) {
               for (int j = -4; j <= 4; j++ ) {
                   if (abs(i) == 4 || abs(j) == 4) {
                       int color = static_cast<int>(image.at<int32_t>(y + i, x + j));
                       if (color == col) {
                           pixel = {(double)x+j-w/2, (double)h/2-(y+i)};
                           return pixel;
                       }
                   }
               }
            }
            return pixel;
        }
    
        bool check_pixel(int x, int y, vector<double> int_point1, vector<double> int_point2, double& error) {
           for (int i = -4; i <= 4; i++ ) {
               for (int j = -4; j <= 4; j++ ) {
                   if (x + j > w || x + j < 0 || y + i > h || y + i < 0) {continue;}
                   int color = static_cast<int>(image.at<int32_t>(y + i, x + j));
                   if (color != 0) {
                       vector<vector<double>> line = img_lines[color];
                       int color = static_cast<int>(image.at<int32_t>(y + i, x + j));
                       vector<double> int_pixel = find_pixel(x + j, y + i, color);
                       vector<double> coord2 = int_pixel;
                       if (coord2[0] < -w) {
                           continue;
                      }
                       vector<double> coord1 = {(double)x+j-w/2, (double)h/2-(y+i)};
                       vector<double> dir1 = direction(line[0], line[line.size()-1]);
                       vector<double> dir2 = direction(int_point1, int_point2);
                       if (check_error(dir1, dir2, error)) {
                           return true;
                       }
                       return false;
                   }
               }
           }
           return false;
        }

        vector<double> direction(vector<double> x, vector<double> y) {
            vector<double> z = sub(x, y);
            double abss = sqrt(z[0]*z[0] + z[1]*z[1]);
            for (int i = 0; i < z.size(); i++) {
                z[i] = z[i] / abss;
            }
            return z;
         }
};




int main() {

   ofstream ply_file("3dlines_test.ply", std::ifstream::out);
   ply_file << "ply\n";
   ply_file << "format ascii 1.0\n";
   ply_file << "comment author: Greg Turk\n";
   ply_file << "comment object: another cube\n";
   ply_file << "element vertex 376\n";
   ply_file << "property float x\n";
   ply_file << "property float y\n";
   ply_file << "property float z\n";
   ply_file << "property uchar red\n";               
   ply_file << "property uchar green\n";
   ply_file << "property uchar blue\n";
   ply_file << "element face 0\n";
   ply_file << "property list uchar int vertex_index\n";
   ply_file << "element edge 0\n";                     
   ply_file << "property int vertex1\n";                
   ply_file << "property int vertex2\n";                 
   ply_file << "property uchar red\n";               
   ply_file << "property uchar green\n";
   ply_file << "property uchar blue\n";
   ply_file << "end_header\n";

   const int cameras_count = 3;
   camera* cameras = new camera[cameras_count];
   cameras[0] = camera(1, "cameras_without_k1.out", "./binary_images/0tt2.png");
   cameras[1] = camera(3, "cameras_without_k1.out", "./test_images/2.png");
   // cameras[2] = camera(2, "cameras_without_k1.out", "./test_images/1.png");
   cameras[2] = camera(4, "cameras_without_k1.out", "./test_images/3.png");
   // cameras[3] = camera(5, "cameras_without_k1.out", "./test_images/4.png");
   // cameras[4] = camera(6, "cameras_without_k1.out", "./test_images/5.png");
   /* cameras[6] = camera(7, "cameras_without_k1.out", "./test_images/6.png");
   cameras[7] = camera(8, "cameras_without_k1.out", "./test_images/7.png");
   cameras[8] = camera(9, "cameras_without_k1.out", "./test_images/8.png");
   cameras[9] = camera(10, "cameras_without_k1.out", "./test_images/9.png");
   cameras[10] = camera(11, "cameras_without_k1.out", "./test_images/10.png"); */
   cout << "cameras" << endl;
   

   int counter = 0;
   for (auto line : cameras[0].get_lines()) {
       if (line.size() < MAX_LINE_SIZE) {
           continue;
       } 
       map<int,int> lines_dict;
       map<int,int> votes;
       vector <map<int, pair<double,int>>> errors(cameras_count - 1);
       vector<map<int,int>> votes_camera(cameras_count - 1);
       vector<double> plane;
       vector<vector<vector<double>>> color_intersections(10000);
       for (int k = 0; k < COUNT_OF_ITERATION; k++) {

           double random_num = rand()*1.0 /  RAND_MAX;
           int index1 = static_cast<int>(line.size() * random_num);
           int index2 = index1 + 5;
           if (index2 >= line.size()) {
               index2 = index1 - 5;
           }

           vector<double> pixel1 = line[index1];
           vector<double> pixel2 = line[index2];
           vector<double> coord1 = cameras[0].change_pixel_to_global_coords(pixel1, 30);
           vector<double> coord2 = cameras[0].change_pixel_to_global_coords(pixel1, 1);
           vector<vector<double>> ep_line;
           cameras[1].find_epipolar_line(coord1, coord2, ep_line);
           vector<vector<vector<double>>> vec(10000);
           int i = 0;
           while (i < ep_line.size()) {
               int x = ep_line[i][0] + cameras[1].get_w() / 2;
               int y = cameras[1].get_h() / 2 - ep_line[i][1];
               int color = static_cast<int>(cameras[1].get_image().at<int32_t>(y, x));
               if (color != 0) {
                   vec[color].push_back({ep_line[i][0], ep_line[i][1]});
               }
               i++;
           }
           vector<vector<double>> points_of_intersect;
           for (int color = 0; color < 10000; color++) {
               if (!vec[color].empty()) {
                   int index = vec[color].size() / 2;
                   points_of_intersect.push_back(vec[color][index]);
               }
           }    
       
           
           for (auto int_point : points_of_intersect) {
               for (int camera_number = 2; camera_number <= cameras_count - 1; camera_number++) {
                   plane = cameras[0].find_plane(pixel1, pixel2);
                   vector<double> pixel_line = cameras[1].find_line(int_point);
                   vector<double> intersect = cameras[0].intersection(plane, pixel_line);
                   vector<double> pixel = {0, 0};
                   try {
                       pixel = cameras[camera_number].change_global_coords_to_pixel(intersect);
                   } catch (string ex) {
                       continue;
                   }
                   int x = static_cast<int>(pixel[0] + cameras[camera_number].get_w()/2);
                   int y = static_cast<int>(cameras[camera_number].get_h()/2 - pixel[1]);
                   int int_x = 0;
                   int int_y = 0;
                   int color = 0;
                   vector<double> int_pixel{-cameras[camera_number].get_w() - 100.0, 0};
                   if (abs(x) < cameras[camera_number].get_w() && abs(y) < cameras[camera_number].get_h()) {
                        int_x = static_cast<int>(int_point[0] + cameras[1].get_w()/2);
                        int_y = static_cast<int>(cameras[1].get_h()/2 - int_point[1]);
                        color = cameras[1].get_image().at<int32_t>(int_y, int_x);
                        vector<vector<double>> int_line = cameras[1].get_lines()[color];
                        int_pixel = cameras[1].find_pixel(int_x, int_y, color);
                        if (int_pixel[0] < -cameras[1].get_w()) {
                            continue;
                        }
                    } else {
                        continue;
                    }
                   
                    auto find_line = cameras[1].get_lines()[color];
                    pixel_line = find_line[0];
                    vector<double> int_point_global1 = cameras[0].intersection(plane, pixel_line);
                    vector<double> int_point1 = cameras[camera_number].change_global_coords_to_pixel(int_point_global1);
                    pixel_line = find_line[find_line.size()-1];
                    vector<double> int_point_global2 = cameras[0].intersection(plane, pixel_line);
                    vector<double> int_point2 = cameras[camera_number].change_global_coords_to_pixel(int_point_global2);
                    double error = 0.0;
                    if (cameras[camera_number].check_pixel(x, y, int_point1, int_point2, error)) {
                        if (lines_dict.count(color) == 0) {
                            lines_dict[color] = 1;
                        } else {
                            lines_dict[color] = lines_dict[color] + 1;
                        }
                        if (errors[camera_number-2].count(color) == 0) {
                            errors[camera_number-2][color].first = 0.0;
                            errors[camera_number-2][color].second = 1;
                        } else {
                            errors[camera_number-2][color].first += error;
                            errors[camera_number-2][color].second += 1;
                        }
                        color_intersections[color].push_back(intersect);
                        votes_camera[camera_number - 2][color] = 1;
                    }
                } // for cameras
              
            } // for int_points
        } // k

        
       /* for (err : errors) {
            cout << "camera " << endl;
            for (auto cur = err.begin(); cur != err.end(); cur++)
	        {
               (*cur).second.first = (*cur).second.first * 1.0 / (*cur).second.second;
               cout << "*** " <<  (*cur).first << " " << (*cur).second.first << endl;
	        }
        } */
        

        for (auto cam_votes : votes_camera) {
            for (auto cur = cam_votes.begin(); cur != cam_votes.end(); cur++)
	        {
               if (votes.count((*cur).first) == 0) {
                   votes[(*cur).first] = 1;
               } else {
                   votes[(*cur).first]++;
               }
	        }
        }

        
        cout << "linnne\n" ;

        for (auto e : errors) {
          cout << "camera" << endl;

        for (auto cur = e.begin(); cur != e.end(); cur++) {
            cout << "err " << (*cur).first << " "  << (*cur).second.first / (*cur).second.second << endl;

        } }

        
        double max = 0.0;
        int key = 0;

        for (auto cur = votes.begin(); cur != votes.end(); cur++)
	    {
            cout << (*cur).first << " " <<  (*cur).second << " " << 1.0*lines_dict[(*cur).first] / (*cur).second << endl;
            if (1.0*lines_dict[(*cur).first] / (*cur).second > max) {
                max = 1.0*lines_dict[(*cur).first] / (*cur).second;
                key = (*cur).first;
            }
            
	    }
 
        cout << "color " << key << endl;
        if (max != 0) {
            for (auto intersect : color_intersections[key]) {
                ply_file << intersect[0] << " " << intersect[1] << " " << intersect[2] << " "  << 255 << " " << 0 << " " << 0 << endl;
                counter++;
            }
        }
     
   }

   cout << counter << endl;

}





