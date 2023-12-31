#include "NURBSCurve.h"

NURBSCurve::NURBSCurve(int _n, int _k, MatrixXd _controlP, VectorXd _knots, bool _isRational)
{
	assert(_k >= 1 && _n >= _k - 1 && _controlP.rows() == _n + 1 && _knots.size() == _n + _k + 1);

	n = _n;
	k = _k;
	knots = _knots;
	controlPw = _controlP;
	isRational = _isRational;

}

bool NURBSCurve::loadNURBS(string name){
	ifstream in(name.c_str());
	if(!in){
		return false;
	}
	in >> isRational;
	in >> n >> k;
	int dimension=3;
	in >> dimension;

	controlPw = MatrixXd(n+1,dimension);
	knots = VectorXd(n+k+1);
	for(int i=0;i < controlPw.rows();i++){
		for(int j=0;j< controlPw.cols();j++){
			in >> controlPw(i,j);
		}
	}
	for(int i=0;i<knots.size();i++){
		in >> knots(i);
	}
	return true;
}

bool NURBSCurve::saveNURBS(string name){
	if(isRational){
		name += ".cptw";
	}else{
		name += ".cpt";
	}
	ofstream out(name.c_str());
	if(!out){
		return false;
	}
	out<< isRational<<endl;
	out<< n <<" "<< k<<endl;
	out<< controlPw.cols()<<endl;
	out<< controlPw <<endl;
	out<< knots.transpose();
	return true;
}

int NURBSCurve::find_ind(double t) const
{
	
	if (t == knots(n + 1)) return n;
	int low = 0;
	int high = n+k;
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

MatrixXd NURBSCurve::eval(double t) const
{
	int L = find_ind(t); 
	assert(L >= k - 1 && L <= n + 1);

	MatrixXd temp = controlPw.block(L - k + 1, 0, k, controlPw.cols());

	for (int r = 1; r <= k - 1; r++)
		for (int i = L - k + 1 + r; i <= L; i++) 
		{
			double factor = (t - knots(i)) / (knots(i + k - r) - knots(i));
			int start = i - (L - k + 1 + r);
			temp.row(start) = (1.0 - factor)*temp.row(start) + factor*temp.row(start + 1);
		}

	return temp.row(0);
}

VectorXd NURBSCurve::parameterize(const MatrixXd & points)
{
	
	int K = points.rows() - 1;
	VectorXd params(points.rows());
	params(0) = 0.0;
	for (int i = 1; i <= K; i++) {
		params(i) = params(i - 1) + (points.row(i) - points.row(i - 1)).norm();
	}
	params = params / params(K);
	params(K) = 1.0;
	return params;
}

double NURBSCurve::basis(int i,int p, double t,const VectorXd &knotvector)
{
	assert(p >= 1);
	if (t<knotvector(i) || t>knotvector(i + p)) {
		return 0.0;
	}
	if (p > 1 && t == knotvector(i + p) && t == knotvector(i + 1)) {
		return 1.0;
	}
	if (p == 1) {
		if (t >= knotvector(i)&& t < knotvector(i+1)) {
			return 1.0;
		}
		else {
			return 0.0;
		}
	}
	
	double a = knotvector(i+p-1) - knotvector(i);
	double b = knotvector(i+p) - knotvector(i+1);
	a = (a == 0.0) ? 0.0 : (t - knotvector(i)) / a;
	b = (b == 0.0) ? 0.0 : (knotvector(i+p) - t) / b;
	return a*basis(i,p-1,t,knotvector) + b*basis(i+1,p-1,t,knotvector);
	
}

void NURBSCurve::interpolate(const MatrixXd &points)
{
	int K = points.rows() - 1;
	VectorXd params(points.rows());
	params(0) = 0.0;
	for (int i = 1; i <= K; i++) {
		params(i) = params(i-1) + (points.row(i) - points.row(i - 1)).norm();
	}
	params = params / params(K);
	params(K) = 1.0;

	knots = VectorXd::Zero(K + 7);
	
	knots.block(3, 0, K + 1, 1) = params;
	knots(K + 6) = 1.0; knots(K + 5) = 1.0; knots(K + 4) = 1.0;
	isRational = false;
	if (points.cols() == 4) {
		isRational = true;
	}
	n = K + 2;
	k = 4;
	controlPw = MatrixXd::Zero(K + 3, points.cols());
	controlPw.row(0) = points.row(0);
	controlPw.row(K + 2) = points.row(K);

	MatrixXd X(K + 1, points.cols()); 
	MatrixXd A = MatrixXd::Zero(K + 1, K + 1);
	
	MatrixXd b = points;

	double delta0 = params(1) - params(0);
	double delta1 = params(2) - params(1);
	double delta2 = params(K) - params(K - 1);
	double delta3 = params(K - 1) - params(K - 2);
	A(0, 0) = (2.0*delta0 + delta1) / (delta0 + delta1);
	A(0, 1) = -delta0 / (delta0 + delta1);
	A(K, K - 1) = -delta2 / (delta2 + delta3);
	A(K, K) = (2.0*delta2 + delta3) / (delta2 + delta3);


	for (int i = 1; i <= K - 1; i++) {
		A(i, i - 1) = basis(i, 4, params(i), knots);
		A(i, i) = basis(i + 1, 4, params(i), knots);
		A(i, i + 1) = basis(i + 2, 4, params(i), knots);

	}
	X = A.inverse()*b;
	controlPw.block(1, 0, K + 1, points.cols()) = X;
	

}

void NURBSCurve::interpolate(const MatrixXd &points, const VectorXd &knotvector)
{
	int K = points.rows() - 1;
	assert(knotvector.size() == points.rows() + 6);
	
	knots = knotvector;
	VectorXd params = knots.block(3, 0, K + 1, 1);
	isRational = false;
	if (points.cols() == 4) {
		isRational = true;
	}
	n = K + 2;
	k = 4;
	controlPw = MatrixXd::Zero(K + 3, points.cols());
	controlPw.row(0) = points.row(0);
	controlPw.row(K + 2) = points.row(K);

	MatrixXd X(K + 1, points.cols());
	MatrixXd A = MatrixXd::Zero(K + 1, K + 1);

	MatrixXd b = points;

	double delta0 = params(1) - params(0);
	double delta1 = params(2) - params(1);
	double delta2 = params(K) - params(K - 1);
	double delta3 = params(K - 1) - params(K - 2);
	A(0, 0) = (2.0*delta0 + delta1) / (delta0 + delta1);
	A(0, 1) = -delta0 / (delta0 + delta1);
	A(K, K - 1) = -delta2 / (delta2 + delta3);
	A(K, K) = (2.0*delta2 + delta3) / (delta2 + delta3);

	for (int i = 1; i <= K - 1; i++) {
		A(i, i - 1) = basis(i, 4, params(i), knots);
		A(i, i) = basis(i + 1, 4, params(i), knots);
		A(i, i + 1) = basis(i + 2, 4, params(i), knots);

	}
	X = A.inverse()*b;
	controlPw.block(1, 0, K + 1, points.cols()) = X;
}

void NURBSCurve::piafit(const MatrixXd &points, int max_iter_num, double eps)
{
	assert(points.rows() > 1 && points.cols() > 0);
	this->n = points.rows() - 1;
	this->k = 4;
	this->isRational = false;
	VectorXd params(points.rows());
	params = parameterize(points);
	cout << "params: " << params.transpose() << endl;
	knots = VectorXd(n + k + 1);
	knots(0) = 0.0; knots(1) = 0.0; knots(2) = 0.0; knots(3) = 0.0;
	knots(n + 4) = 1.0; knots(n + 3) = 1.0; knots(n + 2) = 1.0; knots(n + 1) = 1.0;
	knots.segment(4, n - 3) = params.segment(2, n - 3);
	cout << "knots:" << knots.transpose() << endl;

	controlPw = points;
	double error = 1.0;
	int iter_num = 0;
	MatrixXd delta(points.rows(), points.cols());
	while (error>eps && iter_num<max_iter_num) {
		
		for (int j = 0; j <= n; j++) {
			delta.row(j) = points.row(j) - eval(params(j));
		}
		
		controlPw += delta;
		
		error = delta.rowwise().norm().maxCoeff();
		iter_num++;
		cout << "iter: " << iter_num << ", error: " << error << endl;
	}


}

void NURBSCurve::piafit(const MatrixXd &points, const VectorXd &knotvector,int max_iter_num, double eps)
{
	assert(points.rows() > 1 && points.cols() > 0);
	this->n = points.rows() - 1;
	this->k = 4;
	this->isRational = false;
	VectorXd params(points.rows());
	params = parameterize(points);
	knots = knotvector;


	controlPw = points;
	double error = 1.0;
	int iter_num = 0;
	MatrixXd delta(points.rows(), points.cols());
	while (error>eps && iter_num<max_iter_num) {

		for (int j = 0; j <= n; j++) {
			delta.row(j) = points.row(j) - eval(params(j));
		}

		controlPw += delta;

		error = delta.rowwise().norm().maxCoeff();
		iter_num++;
		cout << "iter: " << iter_num << ", error: " << error << endl;
	}
}

void NURBSCurve::lspiafit(const MatrixXd & points, const int &n_cpts)
{
	assert(points.rows() > 1 && points.cols() > 0);
	this->k = 4;
	const int m = points.rows() - 1;
	this->n = n_cpts - 1;
	const int dimension = points.cols();
	controlPw = MatrixXd(n_cpts, dimension);
	knots = VectorXd(n + k + 1);
	
	VectorXd params = parameterize(points);
	knots(0) = 0.0; knots(1) = 0.0; knots(2) = 0.0; knots(3) = 0.0;
	knots(n + 4) = 1.0; knots(n + 3) = 1.0; knots(n + 2) = 1.0; knots(n + 1) = 1.0;
	double d = 1.0 * (m + 1) / (n - 2);
	for (int j = 1; j <= n - 3; j++) {
		int i = j*d;
		double alpha = 1.0*j*d - i;
		knots(j + 3) = (1.0 - alpha)*params(i - 1) + alpha*params(i);
	}
	cout << "knots: " << knots.transpose() << endl;
	controlPw.row(0) = points.row(0);
	controlPw.row(n) = points.row(m);
	for (int i = 1; i <= n-1; i++) {
		int index = 1.0*(m + 1)*i / n;
		controlPw.row(i) = points.row(index);
	}

	MatrixXd A(m + 1, n + 1);
	for (int i = 0; i <= m; i++) {
		for (int j = 0; j <= n; j++) {
			A(i, j) = basis(j, 4, params(i), knots);
		}
	}
	double mu = 2.0 / (A.transpose()*A).rowwise().sum().maxCoeff();
	const int max_iter_num = 100;
	const double eps = 1e-2;
	double error = 1.0;
	int iter_num = 0;
	
	MatrixXd delta(m + 1, dimension);
	while (error>eps && iter_num<max_iter_num) {
		delta = points - A*controlPw;

		controlPw += mu*A.transpose()*delta;

		error = delta.rowwise().norm().maxCoeff();
		iter_num++;
	}
}

void NURBSCurve::lspiafit(const MatrixXd & points, const int &n_cpts, const VectorXd & knotvector)
{
	assert(points.rows() > 1 && points.cols() > 0);
	this->k = 4;
	const int m = points.rows() - 1;
	this->n = n_cpts - 1;
	const int dimension = points.cols();
	controlPw = MatrixXd(n_cpts, dimension);
	
	VectorXd params = parameterize(points);

	
	knots = knotvector;

	controlPw.row(0) = points.row(0);
	controlPw.row(n) = points.row(m);
	for (int i = 1; i <= n - 1; i++) {
		int index = 1.0*(m + 1)*i / n;
		controlPw.row(i) = points.row(index);
	}
	MatrixXd A(m + 1, n + 1);
	for (int i = 0; i <= m; i++) {
		for (int j = 0; j <= n; j++) {
			A(i, j) = basis(j, 4, params(i), knots);
		}
	}
	double mu = 2.0 / (A.transpose()*A).rowwise().sum().maxCoeff();
	const int max_iter_num = 100;
	const double eps = 1e-2;
	double error = 1.0;
	int iter_num = 0;
	MatrixXd delta(m + 1, dimension);

	while (error>eps && iter_num<max_iter_num) {
		delta = points - A*controlPw;

		controlPw += mu*A.transpose()*delta;

		error = delta.rowwise().norm().maxCoeff();
		iter_num++;
	}
}

bool NURBSCurve::insert(double t)
{
	assert(t >= knots(0) && t <= knots(n + k));
	int L = find_ind(t);
	int start = (L-k+1>=0)?L-k+1:0;
	int end = (L<=n)?L:n;
	MatrixXd new_controlPw(controlPw.rows()+1,controlPw.cols());
	for(int i=0;i<new_controlPw.rows();i++){
		if(i<=start){
			new_controlPw.row(i) = controlPw.row(i);
		}else if(i<=end){
			double factor = (t-knots(i))/(knots(i+k-1)-knots(i));
			new_controlPw.row(i)=factor*controlPw.row(i)+(1.0-factor)*controlPw.row(i-1);
			
		}else{
			new_controlPw.row(i)=controlPw.row(i-1);

		}
	}

	VectorXd new_knots(knots.size()+1);
	for(int i=0;i<new_knots.size();i++){
		if(i<=L){
			new_knots(i)=knots(i);
		}else if(i==L+1){
			new_knots(i)=t;
		}else{
			new_knots(i)=knots(i-1);
		}
	}

	controlPw = new_controlPw;
	knots = new_knots;
	n+=1;

	return true;
}

void NURBSCurve::drawControlPolygon(igl::opengl::glfw::Viewer &viewer){
	
	viewer.data().add_points(controlP, Eigen::RowVector3d(1, 1, 1));
	for (int i = 0; i <= n; i++) {
		viewer.data().add_label(controlP.row(i), std::to_string(i));
	}
	for (int i = 0; i < n; i++){
		viewer.data().add_edges(
			controlP.row(i),
			controlP.row(i + 1),
			Eigen::RowVector3d(1, 1, 1));
	}
}

void NURBSCurve::drawSurface(igl::opengl::glfw::Viewer &viewer, double resolution){
	double left = knots(k - 1);
	double right = knots(n + 1);
	const int num = (right - left) / resolution;
	double new_resolution = (right - left) / num;
	for (int i = 0; i < num; i++)
	{
		Eigen::MatrixXd P1 = eval(left + 1.0*i* new_resolution);
		Eigen::MatrixXd P2 = eval(left + 1.0*(i + 1)*new_resolution);

		if (isRational) {
			viewer.data().add_edges(
				P1.rowwise().hnormalized(),
				P2.rowwise().hnormalized(),
				Eigen::RowVector3d(1, 0, 0));
		}
		else {
			viewer.data().add_edges(P1, P2, Eigen::RowVector3d(1, 0, 0));
		}

	}
}

void NURBSCurve::draw(
	igl::opengl::glfw::Viewer& viewer, 
	bool showpolygon,bool showsurface,
	double resolution)
{
	if(isRational){
		controlP = controlPw.rowwise().hnormalized();
	}else{
		controlP = controlPw;
	}

	if(showpolygon){
		drawControlPolygon(viewer);
	}
	if(showsurface){
		drawSurface(viewer, resolution);
	}
	viewer.core(0).align_camera_center(controlP);
}