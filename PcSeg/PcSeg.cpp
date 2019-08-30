// PcSeg.cpp : 此文件包含 "main" 函数。程序执行将在此处开始并结束。
//相对位置还没考虑
//258的排序方法有点慢
#include "pch.h"
using namespace std;


typedef struct point {
	float x;
	float y;
	float z;

}point ;

typedef struct tile {
	int index;//切块序号
	point tile_xyz[8];//块的八个坐标点

}tile;

//判断该序列的最大包围区域
size_t* getmaxXYZ() {
	pcl::PointCloud<pcl::PointXYZRGB>::Ptr cloud(new pcl::PointCloud<pcl::PointXYZRGB>);

	vector<size_t> x;
	vector<size_t> y;
	vector<size_t> z;
	size_t xyz[6];
	size_t maxX = 0;
	size_t maxY = 0;
	size_t maxZ = 0;
	size_t minX = 10000;
	size_t minY = 10000;
	size_t minZ = 10000;

	for (int j = 1000; j < 1300; j++) {
		if (pcl::io::loadPLYFile<pcl::PointXYZRGB>("loot/ply/loot_vox10_" + to_string(j) + ".ply", *cloud) == -1) //* load the file
		{
			PCL_ERROR("Couldn't read file test_pcd.pcd \n");
		}
		for (size_t i = 0; i < cloud->points.size(); ++i) {
			x.push_back(cloud->points[i].x);
			y.push_back(cloud->points[i].y);
			z.push_back(cloud->points[i].z);
		}
		maxX = (maxX > *max_element(x.begin(), x.end()) ? maxX : *max_element(x.begin(), x.end()));
		maxY = (maxY > *max_element(y.begin(), y.end()) ? maxY : *max_element(y.begin(), y.end()));
		maxZ = (maxZ > *max_element(z.begin(), z.end()) ? maxZ : *max_element(z.begin(), z.end()));
		minX = (minX < *min_element(x.begin(), x.end()) ? minX : *min_element(x.begin(), x.end()));
		minY = (minY < *min_element(y.begin(), y.end()) ? minY : *min_element(y.begin(), y.end()));
		minZ = (minZ < *min_element(z.begin(), z.end()) ? minZ : *min_element(z.begin(), z.end()));
		cout << j << ":" << "xmax: " << maxX << " ymax: " << maxY << " zmax: " << maxZ << " xmin: " << minX << " ymin: " << minY << " zmin: " << minZ << endl;
		x.clear();
		y.clear();
		z.clear();
	}

	cout << maxX << endl << maxY << endl << maxZ << endl << minX << endl << minY << endl << minZ;
	xyz[0] = minX;
	xyz[1] = minY;
	xyz[2] = minZ;
	xyz[3] = maxX;
	xyz[4] = maxY;
	xyz[5] = maxZ;

	return xyz;
}

point* move2x(point *tile_xyz, float x) {
	point temp[8];
	for (int i = 0; i < 8; i++) {
		temp[i].x = tile_xyz[i].x+x;
		temp[i].y = tile_xyz[i].y;
		temp[i].z = tile_xyz[i].z;
	}
	return temp;
}

point* move2y(point *tile_xyz, float y) {

	point temp[8];
	for (int i = 0; i < 8; i++) {
		temp[i].x = tile_xyz[i].x;
		temp[i].y = tile_xyz[i].y+y;
		temp[i].z = tile_xyz[i].z;
	}

	return temp;
}

point* move2z(point *tile_xyz, float z) {

	point temp[8];
	for (int i = 0; i < 8; i++) {
		temp[i].x = tile_xyz[i].x;
		temp[i].y = tile_xyz[i].y;
		temp[i].z = tile_xyz[i].z+z;
	}
	return temp;
}

//自身移动
void movex(tile &tile_xyz, float x) {
	
	for (int i = 0; i < 8; i++) {
		tile_xyz.tile_xyz[i].x += x;
	}
	
}

void movey(tile &tile_xyz, float y) {

	for (int i = 0; i < 8; i++) {
		tile_xyz.tile_xyz[i].y += y;
	}

}

void movez(tile &tile_xyz, float z) {

	for (int i = 0; i < 8; i++) {
		tile_xyz.tile_xyz[i].z += z;
	}

}




//创建对应纵面每层的切块
vector<tile> gettile_seg(int j, const int n, float tile_h, float tile_x, float tile_z) {

	tile tile1;
	tile1.index = j * n*n;
	vector<tile> tiles_seg;

	for (int k = 0; k < 4; k++) {//创建每层的初始块（顶点在(0,y,0)）处
		point point1;
		if (k == 0) {
			point1.x = 0;
			point1.z = 0;
			point1.y = j * tile_h;
			tile1.tile_xyz[k] = point1;
			point1.y = j * tile_h + tile_h;
			tile1.tile_xyz[k + 4] = point1;
		}
		else if (k == 1) {
			point1.x = tile_x;
			point1.z = 0;
			point1.y = j * tile_h;
			tile1.tile_xyz[k] = point1;
			point1.y = j * tile_h + tile_h;
			tile1.tile_xyz[k + 4] = point1;
		}
		else if (k == 2) {
			point1.x = 0;
			point1.z = tile_z;
			point1.y = j * tile_h;
			tile1.tile_xyz[k] = point1;
			point1.y = j * tile_h + tile_h;
			tile1.tile_xyz[k + 4] = point1;
		}
		else if (k == 3) {
			point1.x = tile_x;
			point1.z = tile_z;
			point1.y = j * tile_h;
			tile1.tile_xyz[k] = point1;
			point1.y = j * tile_h + tile_h;
			tile1.tile_xyz[k + 4] = point1;
		}

	}

	tile1.index = j*n*n;
	tiles_seg.push_back(tile1);


	for (int i = 0; i < n; i++) {
		if (i == 0) {
			for (int f = 1; f < n; f++) {
				for (int k = 0; k < 8; k++) {
					tile1.tile_xyz[k] = *(move2x(tiles_seg[f-1].tile_xyz, tile_x) + k);
				}
				tile1.index = f+j*n*n;
				tiles_seg.push_back(tile1);
			}
		}
		else {
			for (int s = 0; s < n; s++) {
				for (int k = 0; k < 8; k++) {
					tile1.tile_xyz[k] = *(move2z(tiles_seg[s+(i-1)*n].tile_xyz, tile_z) + k);
				}
				tile1.index = s + n * i+j*n*n;
				tiles_seg.push_back(tile1);
			}
		}
	}
	return tiles_seg;
}


//获取该序列所有的切块
vector<tile> gettiles(int n, int m, size_t minxyz[3], size_t maxxyz[3], int flag) {//flag=1表示xy平面作为横切面，2表示yz，3表示xz;
	size_t relativexyz[3];
	int tile_number = n * n*m;//块数
	vector<tile> tiles;
	for (int i = 0; i < 3; i++) {
		relativexyz[i] = maxxyz[i] - minxyz[i];
	}
	float relativex = relativexyz[0];
	float relativey = relativexyz[1];
	float relativez = relativexyz[2];

	if (flag == 1) {


	}
	else if (flag == 2) {


	}

	else if (flag == 3) {

		float tile_z = (relativez / n);
		float tile_x = (relativex / n);
		float tile_h = (relativey / m);

		//获取每层的切块
		for (int i = 0; i < m; i++) {

			//获取该层切块
			vector<tile> tile_seg = gettile_seg(i,n, tile_h, tile_x,tile_z);

			for (int k = 0; k < n*n; k++) {
				

				movex(tile_seg[k], minxyz[0]);
				movey(tile_seg[k], minxyz[1]);
				movez(tile_seg[k], minxyz[2]);

				tiles.push_back(tile_seg[k]);
			}
		}




	}


	else {
		//error input
	}

	return tiles;

}



vector<pcl::PointXYZRGB> set_blanktile() {
	pcl::PointXYZRGB point;
	vector<pcl::PointXYZRGB> points;
	point.x = 1;
	point.y = 1;
	point.z = 1;
	point.r = 1;
	point.g = 1;
	point.b = 1;
	for (int i = 0; i < 25; i++) {
		points.push_back(point);

	}
	
	


	return points;
}





//获取所有点对应的切块序号
int gettile_index(vector<pcl::PointXYZRGB, Eigen::aligned_allocator<pcl::PointXYZRGB>>  points, vector<tile> tiles,int n, map <int, vector<pcl::PointXYZRGB>> &index_points) {//n是横截面切分层数
	int ceng = 0;
	size_t tile_number = tiles.size();
	vector<float>::iterator pos;
	

	int y_index;//位置
	int x_index;
	int z_index;
	int index;
	int m = tile_number / n / n;
	size_t maxx = tiles[tile_number - 1].tile_xyz[7].x;
	size_t maxy = tiles[tile_number - 1].tile_xyz[7].y;
	size_t maxz = tiles[tile_number - 1].tile_xyz[7].z;
	size_t minx = tiles[0].tile_xyz[0].x-1;
	size_t miny = tiles[0].tile_xyz[0].y-1;
	size_t minz = tiles[0].tile_xyz[0].z-1;

	size_t tile_h = (maxy-miny) / m;
	size_t tile_x = (maxx-minx) / n;
	size_t tile_z = (maxz - minz) / n;
	size_t x;
	size_t y;
	size_t z;
	vector<float> h;
	vector<float> x_side;
	vector<float> z_side;
	h.push_back(miny);
	x_side.push_back(minx);
	z_side.push_back(minz);
	for (int i = 1; i < m+1; i++) {

		h.push_back(tile_h*i+ miny);

	}
	for (int i = 1; i < n + 1; i++) {
		x_side.push_back(tile_x*i+ minx);
	}

	for (int i = 1; i < n + 1; i++) {
		z_side.push_back(tile_z*i+ minz);
	}


	
	for (size_t i = 0; i < points.size(); i++)
	{
		x = points[i].x;
	    y = points[i].y;
		z = points[i].z;

		//判断属于哪个块
		
		h.push_back(y);
		sort(h.begin(), h.end());
		pos=find(h.begin(), h.end(), y);
		y_index= pos - h.begin();
		h.erase(pos);

		x_side.push_back(x);
		sort(x_side.begin(), x_side.end());
		pos = find(x_side.begin(), x_side.end(), x);
		x_index = pos - x_side.begin();
		x_side.erase(pos);

		z_side.push_back(z);
		sort(z_side.begin(), z_side.end());
		pos = find(z_side.begin(), z_side.end(), z);
		z_index = pos - z_side.begin();
		z_side.erase(pos);
		index = x_index+ n * (z_index - 1) + (n*n)*(y_index - 1);
		

		if (index_points.find(index) !=index_points.end() ) {
			//已存在
			index_points.find(index)->second.push_back(points[i]);
		}
		else {//对应序号不存在 创建一个
			vector<pcl::PointXYZRGB> pointsmap;
			pointsmap.push_back(points[i]);
			index_points.insert(pair<int, vector<pcl::PointXYZRGB>>(index, pointsmap));
		}
	}
	

	return 0;


}

void start(int n, int m ,int j, pcl::PointCloud<pcl::PointXYZRGB>::Ptr cloud) {

	pcl::PointCloud<pcl::PointXYZRGB>::Ptr cloud2(new pcl::PointCloud<pcl::PointXYZRGB>);

	vector<size_t> x;
	vector<size_t> y;
	vector<size_t> z;
	size_t maxX = 0;
	size_t maxY = 0;
	size_t maxZ = 0;
	size_t minX;
	size_t minY;
	size_t minZ;	
	map <int, vector<pcl::PointXYZRGB>> index_points;
	map <int, vector<pcl::PointXYZRGB>>::iterator iter;


		string folderPath = "D:\\dataset\\loot\\loot_"+to_string(n)+"_"+to_string(m)+"\\vox_" + to_string(j);
		if (0 != _access(folderPath.c_str(), 0)) {

			string command;
			command = "mkdir -p " + folderPath;
			system(command.c_str());
		}
		else {
			cout << folderPath << "文件夹已经存在!" << endl;
		}


		for (size_t i = 0; i < cloud->points.size(); ++i) {
			x.push_back(cloud->points[i].x);
			y.push_back(cloud->points[i].y);
			z.push_back(cloud->points[i].z);
		}
		maxX = *max_element(x.begin(), x.end());
		maxY = *max_element(y.begin(), y.end());
		maxZ = *max_element(z.begin(), z.end());
		minX = *min_element(x.begin(), x.end());
		minY = *min_element(y.begin(), y.end());
		minZ = *min_element(z.begin(), z.end());

		size_t minxyz[3] = { minX,minY,minZ };
		size_t maxxyz[3] = { maxX,maxY,maxZ };

		vector<tile> tiles = gettiles(n, m, minxyz, maxxyz, 3);//读取点云文件后获取每个切块边界。


		gettile_index(cloud->points, tiles, n, index_points);//index_points存放每个切块序号对应的点。

		for (int i=1;i<=tiles.size(); i++) {
			cloud2->clear();
			//空切块
			if (index_points.find(i) == index_points.end()) {					
				index_points.insert(pair<int, vector<pcl::PointXYZRGB>>(i, set_blanktile()));	
			}

			size_t pc_number = index_points.find(i)->second.size();//该切块点个数
			cloud2->width = pc_number;
			cloud2->height = 1;
			cloud2->points.resize(cloud2->width * cloud2->height);

			for (size_t k = 0; k < pc_number; k++) {
				cloud2->points[k] = index_points.find(i)->second.at(k);
			}
			pcl::io::savePLYFileASCII(folderPath + "//loot_" + to_string(n) + "-" + to_string(m) + "_" + to_string(j) + "_tile_" + to_string(i) + ".ply", *cloud2);
		}
		index_points.clear();
}

void start2(int n, int m, int j, pcl::PointCloud<pcl::PointXYZRGB>::Ptr cloud) {

	pcl::PointCloud<pcl::PointXYZRGB>::Ptr cloud2(new pcl::PointCloud<pcl::PointXYZRGB>);

	vector<size_t> x;
	vector<size_t> y;
	vector<size_t> z;
	size_t maxX = 0;
	size_t maxY = 0;
	size_t maxZ = 0;
	size_t minX;
	size_t minY;
	size_t minZ;
	map <int, vector<pcl::PointXYZRGB>> index_points;
	map <int, vector<pcl::PointXYZRGB>>::iterator iter;


	string folderPath = "D:\\dataset\\longdress\\longdress_" + to_string(n) + "_" + to_string(m) + "\\vox_" + to_string(j);
	if (0 != _access(folderPath.c_str(), 0)) {

		string command;
		command = "mkdir -p " + folderPath;
		system(command.c_str());
	}
	else {
		cout << folderPath << "文件夹已经存在!" << endl;
	}


	for (size_t i = 0; i < cloud->points.size(); ++i) {
		x.push_back(cloud->points[i].x);
		y.push_back(cloud->points[i].y);
		z.push_back(cloud->points[i].z);
	}
	maxX = *max_element(x.begin(), x.end());
	maxY = *max_element(y.begin(), y.end());
	maxZ = *max_element(z.begin(), z.end());
	minX = *min_element(x.begin(), x.end());
	minY = *min_element(y.begin(), y.end());
	minZ = *min_element(z.begin(), z.end());

	size_t minxyz[3] = { minX,minY,minZ };
	size_t maxxyz[3] = { maxX,maxY,maxZ };

	vector<tile> tiles = gettiles(n, m, minxyz, maxxyz, 3);//读取点云文件后获取每个切块边界。


	gettile_index(cloud->points, tiles, n, index_points);//index_points存放每个切块序号对应的点。

	for (int i = 1; i <= tiles.size(); i++) {
		cloud2->clear();
		//空切块
		if (index_points.find(i) == index_points.end()) {
			index_points.insert(pair<int, vector<pcl::PointXYZRGB>>(i, set_blanktile()));
		}

		size_t pc_number = index_points.find(i)->second.size();//该切块点个数
		cloud2->width = pc_number;
		cloud2->height = 1;
		cloud2->points.resize(cloud2->width * cloud2->height);

		for (size_t k = 0; k < pc_number; k++) {
			cloud2->points[k] = index_points.find(i)->second.at(k);
		}
		pcl::io::savePLYFileASCII(folderPath + "//longdress_" + to_string(n) + "-" + to_string(m) + "_" + to_string(j) + "_tile_" + to_string(i) + ".ply", *cloud2);
	}
	index_points.clear();
}



int main()
{

	pcl::PointCloud<pcl::PointXYZRGB>::Ptr cloud(new pcl::PointCloud<pcl::PointXYZRGB>);


	for (int j = 1000; j < 1300; j++) {

		if (pcl::io::loadPLYFile<pcl::PointXYZRGB>("F:/loot/Ply/loot_vox10_" + to_string(j) + ".ply", *cloud) == -1) //* load the file
		{

			PCL_ERROR("Couldn't read file\n");

		}
		start(2, 2, j, cloud);
		start(2, 3, j, cloud);
		start(3, 2, j, cloud);
		start(3, 3, j, cloud);
		start(4, 2, j, cloud);
		start(4, 3, j, cloud);
	}

	for (int j = 1051; j < 1350; j++) {

		if (pcl::io::loadPLYFile<pcl::PointXYZRGB>("F://longdress/Ply/longdress_vox10_" + to_string(j) + ".ply", *cloud) == -1) //* load the file
		{

			PCL_ERROR("Couldn't read file\n");

		}

		start2(2, 2, j, cloud);
		start2(2, 3, j, cloud);
		start2(3, 2, j, cloud);
		start2(3, 3, j, cloud);
		start2(4, 2, j, cloud);
		start2(4, 3, j, cloud);

	}







}