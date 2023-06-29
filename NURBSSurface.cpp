#include "NURBSSurface.h"

NURBSSurface::NURBSSurface(
	VectorXi _order, 
	vector<MatrixXd> _controlP, 
	VectorXd _uknots, 
	VectorXd _vknots,
	bool _isRational)
{
	u_order = _order(0);
	v_order = _order(1);
	v_num = _controlP.size() - 1;
	assert(u_order >= 1 && v_order >= 1 && v_num >= 1);
	u_num = _controlP[0].rows() - 1;
	assert(_uknots.size() == u_num + u_order + 1 && _vknots.size() == v_num + v_order + 1);
	assert(u_num >= u_order - 1 && v_num >= v_order - 1);

	dimension = _controlP[0].cols();
	uknots = _uknots;
	vknots = _vknots;
	controlPw = _controlP;
	isRational = _isRational;

	if (isRational) { assert(dimension == 4); }
	else { assert(dimension == 3); }

}

bool NURBSSurface::loadNURBS(string name){
	ifstream in(name);
	if(!in){
		return false;
	}
	char sep;
	dimension=3;
	in>> isRational;
	in>> u_num >> v_num;
	in>> u_order >> v_order;
	in>> dimension;
	controlPw = vector<MatrixXd>(v_num+1);
	for(int i=0;i<controlPw.size();i++){
		controlPw[i] = MatrixXd(u_num+1,dimension);
	}
	for(int i=0;i<controlPw.size();i++){
		for(int j=0;j<=u_num;j++){
			for(int k=0;k<dimension;k++){
				in>> controlPw[i](j,k);
			}
			in>>sep;
		}
	}

	uknots = VectorXd(u_num+u_order+1);
	vknots = VectorXd(v_num+v_order+1);

	for(int i=0;i<uknots.size();i++){
		in>>uknots(i);
	}
	
	for(int i=0;i<vknots.size();i++){
		in>>vknots(i);
	}
	return true;
}

bool NURBSSurface::saveNURBS(string name){
	if(controlPw.size()==0){
		cout<< "nothing to save!"<<endl;
		return false;
	}
	if(isRational){
		name+=".cptw";
	}else{
		name+=".cpt";
	}
	ofstream out(name);
	if(!out){
		return false;
	}
	IOFormat outputFmt(4, 0, " ", " ", "", ",");	
	out<< isRational<<endl;
	out<< u_num << " " <<v_num<<endl;
	out<< u_order << " "<< v_order<<endl;
	out<< controlPw[0].cols()<<endl;
	
	for(int i=0;i<controlPw.size();i++){
		for(int j=0;j<=u_num;j++){
			out<<controlPw[i].row(j).format(outputFmt);
		}
		out<<endl;
	}
	out<<uknots.transpose()<<endl;
	out<<vknots.transpose();
	return true;
}
	
int NURBSSurface::find_ind(double t, int k, int n, const VectorXd& knots)
{
	if (t == knots(n + 1)) return n;
	int low = 0;
	int high = n + k;
	assert(t >= knots(low) && t < knots(high));

	int mid = (low + high) / 2;
	while (t < knots(mid) || t >= knots(mid + 1))
	{
		if (t < knots(mid)) high = mid;
		else low = mid;
		mid = (low + high) / 2;
	}
	return mid;
}

MatrixXd NURBSSurface::eval(double u, double v)
{
	MatrixXd v_controlP(v_num + 1, dimension); 
	for (int i = 0; i <= v_num; i++)
	{
		v_controlP.row(i) = eval(u, controlPw[i], uknots);
	}
	return eval(v, v_controlP, vknots); 
}

MatrixXd NURBSSurface::eval(
	double t, 
	const MatrixXd &_controlP, 
	const VectorXd &knots)
{
	
	int n = _controlP.rows() - 1;
	int k = knots.size() - _controlP.rows();
	assert(t>=knots(k-1) && t<=knots(n+1));
	int L = find_ind(t, k, n, knots); 

	MatrixXd temp = _controlP.block(L - k + 1, 0, k, _controlP.cols());

	for (int r = 1; r <= k - 1; r++)
		for (int i = L - k + 1 + r; i <= L; i++) 
		{
			double factor = (t - knots(i)) / (knots(i + k - r) - knots(i));
			int start = i - (L - k + 1 + r);
			temp.row(start) = (1.0 - factor)*temp.row(start) + factor*temp.row(start + 1);
		}
	MatrixXd curvePoint = temp.row(0);
	return curvePoint;
}

bool NURBSSurface::insert(double s, char dir){
	assert(dir=='u' || dir=='v');
	if(dir=='u'){
		for(int i=0;i<controlPw.size();i++){
			NURBSCurve nurbs(u_num,u_order,controlPw[i],uknots,isRational);
			nurbs.insert(s);
			controlPw[i]=nurbs.controlPw;
			if (i == controlPw.size() - 1) {
				uknots = nurbs.knots;
				u_num = nurbs.n;
			}
		}
		return true;
	}else if(dir=='v'){
		vector<MatrixXd> new_controlPw(controlPw.size()+1);
		for (int i = 0; i < new_controlPw.size(); i++) {
			new_controlPw[i] = MatrixXd(u_num + 1, controlPw[0].cols());
		}
		MatrixXd v_controlPw(v_num+1,controlPw[0].cols());
		for(int i=0;i<=u_num;i++){
			for(int j=0;j<=v_num;j++){
				v_controlPw.row(j) = controlPw[j].row(i);
			}
			NURBSCurve nurbs(v_num,v_order,v_controlPw,vknots,isRational);
			nurbs.insert(s);
			
			for(int k=0;k<new_controlPw.size();k++){
				new_controlPw[k].row(i) = nurbs.controlPw.row(k);
			}
			if (i == u_num) {
				vknots = nurbs.knots;
				v_num = nurbs.n;
			}
		}
		controlPw = new_controlPw;
		return true;

	}else{
		cout<< "please input dir as u or v!"<<endl;
		return false;
	}
}

bool NURBSSurface::insert(double s, double t){
	if(insert(s,'u') && insert(t,'v')){
		return true;
	}
	return false;
}

void NURBSSurface::drawControlPolygon(igl::opengl::glfw::Viewer &viewer){
	for (int i = 0; i < controlP.size(); i++)
	{
		viewer.data().add_points(
			controlP[i],
			Eigen::RowVector3d(1, 1, 1));
	}
	for (int j = 0; j <= v_num; j++)
		for (int i = 0; i <= u_num; i++)
		{
			if (i != u_num) {
				viewer.data().add_edges(
					controlP[j].row(i),
					controlP[j].row(i + 1),
					Eigen::RowVector3d(1, 1, 1));
			}
			if (j != v_num) {
				viewer.data().add_edges(
					controlP[j].row(i),
					controlP[j + 1].row(i),
					Eigen::RowVector3d(1, 1, 1));
			}
		}
}

void NURBSSurface::drawSurface(igl::opengl::glfw::Viewer &viewer, double resolution){
	double u_low = uknots(u_order - 1);
	double u_high = uknots(u_num + 1);
	const int uspan = (u_high - u_low) / resolution;
	double u_resolution = (u_high - u_low) / uspan;

	double v_low = vknots(v_order - 1);
	double v_high = vknots(v_num + 1);
	const int vspan = (v_high - v_low) / resolution;
	double v_resolution = (v_high - v_low) / vspan;

	mesh_V = MatrixXd((uspan + 1)*(vspan + 1), 3);
	mesh_F = MatrixXi(2 * uspan*vspan, 3);
	for (int j = 0; j <= vspan; j++)
		for (int i = 0; i <= uspan; i++)
		{
			RowVectorXd curvePoint = eval(u_low + i*u_resolution, v_low + j*v_resolution).row(0);
			if (isRational) { mesh_V.row(j*(uspan + 1) + i) = curvePoint.hnormalized(); }
			else { mesh_V.row(j*(uspan + 1) + i) = curvePoint; }
		}

	for (int j = 0; j<vspan; j++)
		for (int i = 0; i < uspan; i++)
		{
			int V_index = j*(uspan + 1) + i;
			int F_index = 2 * j*uspan + 2 * i;
			mesh_F.row(F_index) << V_index, V_index + 1, V_index + uspan + 1;
			mesh_F.row(F_index + 1) << V_index + uspan + 1, V_index + 1, V_index + uspan + 2;
		}

	viewer.data().set_mesh(mesh_V, mesh_F);
}


void NURBSSurface::draw(
	igl::opengl::glfw::Viewer &viewer, 
	bool showpolygon,bool showsurface,
	double resolution){

	if(controlP.size()!=controlPw.size()){
		controlP = vector<MatrixXd>(controlPw.size());
	}
	if(isRational){
		for(int i=0;i<controlP.size();i++){
			controlP[i] = controlPw[i].rowwise().hnormalized();
		}
	}else{
		controlP = controlPw;
	}
	
	if(showpolygon){
		drawControlPolygon(viewer);
		viewer.core(0).align_camera_center(controlP[controlP.size()/2]);
	}
	if(showsurface){
		drawSurface(viewer,resolution);
		viewer.core(0).align_camera_center(mesh_V,mesh_F);
	}

}

void NURBSSurface::skinning(const vector<NURBSCurve> &curves,igl::opengl::glfw::Viewer &viewer){
	assert(curves.size()>1);
	isRational = curves[0].isRational;
	vector<map<double,int>> curve_knots(curves.size());
	for(int i=0;i<curves.size();i++){
		for(int j=0;j<curves[i].knots.size();j++){
			if(curve_knots[i].find(curves[i].knots(j))==curve_knots[i].end()){
				curve_knots[i][curves[i].knots(j)] = 1;
			}else{
				curve_knots[i][curves[i].knots(j)] += 1;
			}
		}
	}
	map<double, int> merge_knots;
	for (int i = 0; i < curve_knots.size(); i++) {
		for (auto it = curve_knots[i].begin(); it != curve_knots[i].end(); it++) {
			if (merge_knots.find(it->first) == merge_knots.end()) {
				merge_knots[it->first] = it->second;
			}
			else if (merge_knots[it->first] < it->second) {
				merge_knots[it->first] = it->second;
			}
		}
	}
	cout << "merge_knots:\n" << endl;
	for (auto it = merge_knots.begin(); it != merge_knots.end(); it++) {
		cout << it->first << ", " << it->second << endl;
	}

	vector<NURBSCurve> new_curves = curves;

	for (int i = 0; i < curve_knots.size(); i++) {
		for (auto it = merge_knots.begin(); it != merge_knots.end(); it++) {
			int insert_num = 0;
			if (curve_knots[i].find(it->first) == curve_knots[i].end()) {
				insert_num = it->second;
			}
			else {
				insert_num = it->second - curve_knots[i][it->first];
			}
			for (int j = 0; j < insert_num; j++) {
				new_curves[i].insert(it->first);
			}
		}
	}

	
	for (int i = 0; i < new_curves.size(); i++) {
		cout << "curve " << i << ": " << new_curves[i].knots.transpose() << endl;
		new_curves[i].draw(viewer);
	}
	

	vknots = new_curves[0].knots;
	VectorXd curves_param = VectorXd::Zero(new_curves.size());
	v_num = new_curves[0].n;
	v_order = new_curves[0].k;
	dimension = new_curves[0].controlPw.cols();
	u_order = 4;
	int curve_num = curves.size() - 1;
	u_num = curve_num + 2;

	for (int i = 0; i <= v_num; i++) {
		MatrixXd u_cpts(curve_num + 1, dimension);
		for (int j = 0; j <= curve_num; j++) {
			u_cpts.row(j) = new_curves[j].controlPw.row(i);
		}
		curves_param += NURBSCurve::parameterize(u_cpts);
	}
	curves_param /= (v_num + 1);
	uknots = VectorXd::Zero(curve_num + 7);

	uknots.block(3, 0, curve_num + 1, 1) = curves_param;
	uknots(curve_num + 6) = 1.0; uknots(curve_num + 5) = 1.0; uknots(curve_num + 4) = 1.0;

	cout << "uknots: " << uknots.transpose() << endl;
	controlPw = vector<MatrixXd>(v_num + 1);
	for (int i = 0; i < controlPw.size(); i++) {
		controlPw[i] = MatrixXd(u_num + 1, dimension);
	}

	for (int i = 0; i <= v_num; i++) {
		MatrixXd u_ctps = MatrixXd::Zero(curve_num + 1, dimension);

		for (int j = 0; j <= curve_num; j++) {
			u_ctps.row(j) = new_curves[j].controlPw.row(i);
			
		}
		NURBSCurve nurbs;
		
		nurbs.interpolate(u_ctps, uknots);
		nurbs.draw(viewer, false, true);
		controlPw[i] = nurbs.controlPw;

	}
	cout << "skinning finished!" << endl;
}
