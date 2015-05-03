#include <iostream>
#include <vector>
#include <map>
#include <algorithm>
#include <set>
#include <queue>
#include <stack>
#include <string>
#include <sstream>
#include <cstdio>
#include <cstring>
#include <cmath>
#include <cassert>
#include <time.h>     

using namespace std;

#define ri(X) scanf("%d", &(X))
#define pi(X) printf("%d", (X))
#define mp(X,Y) make_pair(X,Y)
#define pb(X) push_back(X)
#define lint long long
#define pii pair<int,int>
#define inf 1e9
#define linf 1e18
#define X first
#define Y second
#define all(X) (X).begin(),(X).end()
#define uni(X) X.erase(unique(X.begin(), X.end()), X.end());

double INF = 1e100;
double EPS = 1e-12;






string itos(int x){
  stringstream ss;
  ss << x;
  return ss.str();
}


const double pi = acos(-1);
class PT{
  public:
  double x, y;
  int i;
  PT(int xx = 0, int yy = 0, int ii = 0){
    x = xx;
    y = yy;
    i = ii;
  }
  bool operator <(const PT &p)  const { return x < p.x || (x == p.x && y < p.y); }
  bool operator ==(const PT &p)  const { return x == p.x && y == p.y; }
  PT operator + (const PT &p)  const { return PT(x+p.x, y+p.y); }
  PT operator - (const PT &p)  const { return PT(x-p.x, y-p.y); }
  PT operator * (double c)     const { return PT(x*c,   y*c  ); }
  PT operator / (double c)     const { return PT(x/c,   y/c  ); }
};
double dot(PT p, PT q)     { return p.x*q.x+p.y*q.y; }
double dist2(PT p, PT q)   { return dot(p-q,p-q); }
double dist(PT p, PT q)   { return sqrt(dist2(p,q)); }
double cross(PT p, PT q)   { return p.x*q.y-p.y*q.x; }

double angle0(PT p){
  return atan2(p.x,p.y);
}

ostream &operator<<(ostream &os, const PT &p) {
  os << "(" << p.x << "," << p.y << ")"; 
}
vector<PT> hull(vector<PT> P, bool b = false){
  assert(P.size() >= 3);
	int n = P.size(), k = 0; 
	vector<PT> H(2*n); 
	// Sort PTs lexicographically 
	sort(P.begin(), P.end()); 
	// Build lower hull 
	for (int i = 0; i < n; i++) {
	  if(b == 0){
		  while (k >= 2 && cross(H[k-1]-H[k-2], P[i]-H[k-2]) < 0) k--; 
		}else{
		  while (k >= 2 && cross(H[k-1]-H[k-2], P[i]-H[k-2]) > 0) k--; 
		}
		H[k++] = P[i]; 
	}
	H.resize(k);
	return H; 
}


vector<PT> convex_hull(vector<PT> P){ 
	int n = P.size(), k = 0; 
	vector<PT> H(2*n); 
	// Sort PTs lexicographically 
	sort(P.begin(), P.end()); 
	// Build lower hull 
	for (int i = 0; i < n; i++) { 
		while (k >= 2 && cross(H[k-1]-H[k-2], P[i]-H[k-2]) < 0) k--; 
		H[k++] = P[i]; 
	} 
	// Build upper hull 
	for (int i = n-2, t = k+1; i >= 0; i--) { 
		while (k >= t && cross(H[k-1]-H[k-2], P[i]-H[k-2]) < 0) k--; 
		H[k++] = P[i]; 
	} 
	H.resize(k);
	H.pop_back();
	return H; 
}

// O(ni)
double ComputeArea(const vector<PT> &p) {
  double area = 0;
  for(int i = 0; i < p.size(); i++) {
    int j = (i+1) % p.size();
    area += p[i].x*p[j].y - p[j].x*p[i].y;
  }
  return fabs(area / 2.0);
}
// O(1)
bool LinesParallel(PT a, PT b, PT c, PT d) { 
  return fabs(cross(b-a, c-d)) < EPS; 
}

// O(1)
bool LinesCollinear(PT a, PT b, PT c, PT d) { 
  return LinesParallel(a, b, c, d)
      && fabs(cross(a-b, a-c)) < EPS
      && fabs(cross(c-d, c-a)) < EPS; 
}

// O(1)
bool SegmentsIntersect(PT a, PT b, PT c, PT d) {
  if (LinesCollinear(a, b, c, d)) {
    if (dist2(a, c) < EPS || dist2(a, d) < EPS ||
      dist2(b, c) < EPS || dist2(b, d) < EPS) return true;
    if (dot(c-a, c-b) > 0 && dot(d-a, d-b) > 0 && dot(c-b, d-b) > 0)
      return false;
    return true;
  }
  if (cross(d-a, b-a) * cross(c-a, b-a) > 0) return false;
  if (cross(a-c, d-c) * cross(b-c, d-c) > 0) return false;
  return true;
}

// O(ni)
bool hasArea(const vector<PT> &p){
  bool allcol = 1;
  for(int i = 0; i < p.size()-1 && allcol; i++){
    if(!LinesParallel(p[i], p[i+1], p[i+1], p[(i+2)%p.size()])){
      allcol = 0;
    }
  }
  if(allcol) return 0;
  return 1;
}
// O(ni*ni)
bool IsSimple(const vector<PT> &p) {
  return true;
  
  for (int i = 0; i < p.size(); i++) {
    for (int k = i+1; k < p.size(); k++) {
      int j = (i+1) % p.size();
      int l = (k+1) % p.size();
      if (i == l || j == k) continue;
      if (SegmentsIntersect(p[i], p[j], p[k], p[l])){
        cerr << p[i] << p[j] << p[k] << p[l] << endl;
        return false;
      }
    }
  }
  return true;
}
vector<PT> ve;
vector<string> res;
int np;
int n;

vector<int> reorder(PT center, vector<int> vi){
  int ki = 0;
  for(int i = 0; i < vi.size(); i++){
    if(ve[vi[i]] == center){
      ki = i;
    }
  }
  if(ki == 0) return vi;
  rotate(vi.begin(), vi.begin()+ki, vi.end());
  
  assert(ve[vi[0]] == center);
  
  return vi;
}


double vicost(vector<int> vi){
  vector<PT> vq;
  double cost;
  for(int i = 0; i < vi.size(); i++){
    vq.pb(ve[vi[i]]);
  }
  if(!hasArea(vq) || !IsSimple(vq)){
    //cerr << "SHiT" << endl;
    cost = inf;
    return cost;
  }
  cost = ComputeArea(vq);
  return cost;
}


vector<int> getpol1(vector<int>, int, int);
vector<int> getpol(vector<int> vi);
vector<int> getpol2(vector<int> , bool);
vector<int> getpol3(vector<int> vi);


class poly{
public:
  vector<int> vi;
  double cost;
  poly(vector<int> vv){
    cost = -1;
    vi = vv;
  }
  poly(){
    cost = inf;
  }  
  // O(C(Ri))
  double getCost(){
    if(cost >= 0){
      return cost;
    }
    return cost = vicost(vi);
  }
  // O(npi)
  string tostring(){
    string rp;
    for(int i = 0; i < vi.size(); i++){
      rp+=itos(vi[i])+ " ";
    }
    return rp;
  }
  
  bool improvepol(){
    vector<pair<double,vector<int> > > vs;
    vector<int> p;
    double c;
    
    
    p = getpol3(vi);
    c = vicost(p);
    vs.pb(mp(c,p));
    
    /*p = getpol2(vi,0);
    c = vicost(p);
    vs.pb(mp(c,p));
    
    for(int i = 0; i < 2; i++){
      p = getpol2(vi,1);
      c = vicost(p);
      vs.pb(mp(c,p));
    }*/
    
    for(int i = 0; i < 7; i++){
      p = getpol(vi);
      c = vicost(p);
      vs.pb(mp(c,p));
    }
    
    int mi = 0;
    for(int i = 0; i < vs.size(); i++){
      if(vs[i].X < vs[mi].X){
        mi = i;
      }
    } 
    
    if(getCost() <= vs[mi].X){
      return false;
    }
    
    cost = vs[mi].X;
    vi = vs[mi].Y;
    return true;
  }
  
  
};




// O(ni lg ni + C(R));
vector<poly> split2(vector<int> vi, int ng){
  
  assert(vi.size()  >=  6);
  
  vector<int> vr;
  vector<bool> visit;
  vector<int> match;
  for(int i = 0;i < vi.size(); i++){
    vr.pb(i);
    visit.pb(0);
    match.pb(0);
  }
  random_shuffle(vr.begin(), vr.end());
  
  
  vector<vector<int> > vs;
  
  vector<poly> res;
  
  vector<PT> vp;
  for(int i = 0; i < min((int)vi.size()/3,ng); i++){
    vector<int> vini;
    vini.pb(vi[vr[i]]);
    vs.pb(vini);
    vp.pb(ve[vi[vr[i]]]);
    visit[vr[i]] = 1;
    match[vr[i]] = i;
  }
  
  
  for(int i = 0; i < vi.size(); i++){
    if(visit[i]) continue;
    for(int j = 0; j < vs.size(); j++){
      PT mnew = vp[j]/vs[j].size();
      PT mold = vp[match[i]]/vs[match[i]].size();
      if(dist2(ve[vi[i]], mnew) < dist2(ve[vi[i]], mold)){
        match[i] = j;
      }
    }
  }
  
  for(int k = 0; k < vi.size() - vs.size(); k++){
    
    int mi = 0;
    double mid = inf;
    for(int i = 0; i < vi.size(); i++){
      if(visit[i]) continue;
      
      PT mean = vp[match[i]]/vs[match[i]].size();
      
      double d = dist2(ve[vi[i]], mean);
      if(d < mid){
        d = mid;
        mi = i;
      }
    }
    
    visit[mi] = 1;
    vs[match[mi]].pb(vi[mi]);
    vp[match[mi]] = vp[match[mi]]+ ve[vi[mi]];
    
    for(int i = 0; i < vi.size(); i++){
      if(visit[i]) continue;
      PT mnew = vp[match[mi]]/vs[match[mi]].size();
      PT mold = vp[match[i]]/vs[match[i]].size();
      if(dist2(ve[vi[i]], mnew) < dist2(ve[vi[i]], mold)){
        match[i] = match[mi];
        //cerr <<"hola" << endl;
      }
      for(int j = 0; j < vs.size(); j++){
        PT mnew = vp[j]/vs[j].size();
        PT mold = vp[match[i]]/vs[match[i]].size();
        if(dist2(ve[vi[i]], mnew) < dist2(ve[vi[i]], mold)){
          match[i] = j;
        }
        //assert(dist2(ve[vi[i]], mnew) >= dist2(ve[vi[i]], mold));
      }
    }
  }
  
  for(int i = 0; i < vs.size(); i++){
    vector<PT> vq;
    for(int j = 0; j < vs[i].size(); j++){
      vq.pb(ve[vs[i][j]]);
    }
    if(vs.size() < 3 || !hasArea(vq)){
      return res;
    }
  }
  
  for(int i = 0; i < vs.size(); i++){
    res.pb(getpol(vs[i]));
  }
  
  return res;
  
}


// O(ni lg ni + C(R));
pair<poly, poly> split(vector<int> vi, int w){
  
  assert(vi.size()  >=  6);
  //cerr <<"Hola:  "<< k << " " << vp.size() <<" " << endl;
  int rot = rand()%2;
  vector<pair<pair<int,int>,int> > vs;
  for(int i = 0; i < vi.size(); i++){
    if(rot == 0){
      vs.pb(mp(mp(ve[vi[i]].x,ve[vi[i]].y),vi[i]));
    }else{
      vs.pb(mp(mp(ve[vi[i]].y,ve[vi[i]].x),vi[i]));
    }
  }
  sort(all(vs));
  vector<int> rp[2];
  for(int i = 0; i < vs.size(); i++){
    rp[i< w?0:1].pb(vs[i].Y);//
  }
  
  return mp(getpol(rp[0]), getpol(rp[1]));
}

class solution{
  public:
  vector<poly> vp;
  double cost;
  solution(){
    cost = inf;
  }
  solution(vector<poly> vv){
    vp = vv;
    cost = -1;
  }
  
  // O(C(R)) 
  double getCost(){
    if(cost >= 0){
      return cost;
    }
    cost = 0;
    for(int i = 0; i < vp.size(); i++){
      cost+=vp[i].getCost();
    }
    return cost;
  }
  
  // O(np)
  vector<string> tovs(){
    vector<string> vs;
    for(int i = 0; i < vp.size(); i++){
      vs.pb(vp[i].tostring());
    }
    return vs;
  }
  
  
  
  // O(ni lg ni + c(R))
  void splitr(int k){
    if(vp.size() >= n) return;
    assert(0 <= k && k < vp.size());
    int siz = vp[k].vi.size();
    if(siz <= 8) return;
    
    
    
    pair<poly, poly> vn = split(vp[k].vi, 3+(rand()%(siz-6)));
    
    vp[k] = vn.X;
    vp.pb(vn.Y);
    
    cost = -1;
  }
  
  // O(ni lg ni + C(R))
  void star(int k){
    assert(0 <= k && k < vp.size());
    vp[k] = getpol(vp[k].vi);
    cost = -1;
  }
  
  // O(ni lg ni + C(R))
  void semistar(int k){
    assert(0 <= k && k < vp.size());
    vp[k] = getpol2(vp[k].vi, -1);
    cost = -1;
  }
  
  // O(ni lg ni + C(R))
  void implosion(int k){
    assert(0 <= k && k < vp.size());
    vp[k] = getpol3(vp[k].vi);
    cost = -1;
  }
  
  
  bool improvepol(int k){
    
    bool res = vp[k].improvepol();
    if(res == 1) cost = -1;
    return res;
  }
  
  bool divide(int k) {
    if(vp.size() >= n) return false;
    assert(0 <= k && k < vp.size());
    
    int siz = vp[k].vi.size();
    if(siz <= 6) return false;
    bool changed = false;
    pair<poly,poly > vn, best;
    int div;
    
    
    
    for(int i = 0; i < 7; i++){
      if(i == 0){
        div = siz/2;
      }else{
        div = 3+(rand()%(siz-6));
      }
      
      vn = split(vp[k].vi, div);
      vn.X.improvepol();
      vn.Y.improvepol();
      if(vp[k].getCost() <= vn.X.getCost()+vn.Y.getCost()) continue;
      
      best = vn;
      changed = true;
    }
    
    if(!changed) return false;
    
    vp[k] = best.X;
    vp.pb(best.Y);
    cost = -1;
    
    return true;
  }
  bool divide2() {
    
    assert(vp.size() == 1);
    int siz = vp[0].vi.size();
    if(siz <= 6) return false;
    
    
    vector<poly> res;
    double sum = 0;
    res = split2(vp[0].vi, n);
    
    if(res.size() <= 1) return 0;
    
    
    
    for(int i = 0; i < res.size(); i++){
      res[i].improvepol();
      sum+= res[i].getCost();
    }
    
    if(vp[0].getCost() <= sum) return 0;
    
    vp.pop_back();
    for(int i = 0; i < res.size(); i++){
      vp.pb(res[i]);
    }
    cost = sum;
    return true;
    
  }
  
  
};

vector<int> fixend(vector<int> vi){
  vector<pair<double, int> > vc;
  int si = vi.size()-1;
  for(int k = vi.size()-1; k >= 0; k--){
    if(LinesParallel(ve[vi[0]], ve[vi[vi.size()-1]], ve[vi[(k+1)%vi.size()]], ve[vi[k]])){
      double di = dist2(ve[vi[0]], ve[vi[k]]);
      vc.pb(mp(-di,vi[k]));
      si = k;
    }else{
      break;
    }
  }
  sort(all(vc));
  for(int i = 0; i < vc.size(); i++){
    vi[si+i] = vc[i].Y;
  }
  return vi;
}

// O(ni lg ni)
vector<int> getpol1(vector<int> vi, int s=0){
  assert(vi.size() >= 2);
  if(vi.size() == 3) return vi;
  PT ic = ve[vi[s]];
  
  vector<pair<pair<double,double>, int> > vs;
  for(int i = 0; i < vi.size(); i++){
    if(vi[i] == vi[s]){
      continue;
    }
    double ang = angle0(ve[vi[i]]-ic);
    double di = dist2(ic, ve[vi[i]]);
    vs.pb(mp(mp(ang,di),vi[i]));
  }
  sort(all(vs));
  vi[0] = vi[s];
  
  int ki = 0;
  for(int i = 0; i < vs.size()-1; i++){
    if(fabs(vs[i].X.X-vs[i+1].X.X) >= pi){
      ki = i+1;
    }
  }
  for(int i = 0; i < vs.size(); i++){
    vi[i+1] = vs[(i+ki)%vs.size()].Y;
  }
  
  vi = fixend(vi);
  return (vi);
}





// O(ni lg ni + C(R));
vector<int> getpol2(vector<int> vi, bool q = 0){
  if(vi.size() <= 8) return getpol(vi);
  int siz = vi.size();
  int w;
  if(q == 0){
    w =vi.size()/2;
  }else{
    w = 3+(rand()%(siz-6));
  }
  //cerr <<"Hola:  "<< k << " " << vp.size() <<" " << endl;
  int rot = 1;//rand()%2;
  vector<pair<pair<int,int>,int> > vs;
  for(int i = 0; i < vi.size(); i++){
    if(rot == 0){
      vs.pb(mp(mp(ve[vi[i]].x,ve[vi[i]].y),vi[i]));
    }else{
      vs.pb(mp(mp(ve[vi[i]].y,ve[vi[i]].x),vi[i]));
    }
  }
  sort(all(vs));
  vector<int> rp[2];
  vector<PT> vh[2];
  for(int i = 0; i < vs.size(); i++){
    rp[i< w?0:1].pb(vs[i].Y);//
    vh[i< w?0:1].pb(ve[vs[i].Y]);//
  }
  //cout << w << endl;
  vh[0] = hull(vh[0],1);
  vh[1] = hull(vh[1],0);
  
  
  rp[0] = reorder(vh[0][1], rp[0]);
  rp[1] = reorder(vh[1][1], rp[1]);
  
  

  rp[0] = getpol1(rp[0],0);
  rp[1] = getpol1(rp[1],0);
  
  reverse(all(rp[1]));
  
  rp[0] = reorder(vh[0][0], rp[0]);
  rp[1] = reorder(vh[1][0], rp[1]);
  
  
  if(
     SegmentsIntersect(ve[rp[0][0]],ve[rp[1][0]], ve[rp[0][1]],ve[rp[1][1]]) || 
     SegmentsIntersect(ve[rp[0][0]],ve[rp[1][0]], ve[rp[0][1]],ve[rp[0][2]]) ||
     SegmentsIntersect(ve[rp[0][0]],ve[rp[1][0]], ve[rp[1][1]],ve[rp[1][2]]) ||
     SegmentsIntersect(ve[rp[0][1]],ve[rp[1][1]], ve[rp[0][0]],ve[rp[0][2]]) ||
     SegmentsIntersect(ve[rp[0][1]],ve[rp[1][1]], ve[rp[1][0]],ve[rp[1][2]])
     ){
    return getpol1(vi);
  }
  int k = 0;
  for(int i = 1; i < rp[0].size(); i++){
    vi[k++] = rp[0][i];
  }
  vi[k++] = rp[0][0];
  vi[k++] = rp[1][0];
  for(int i = rp[1].size()-1; i >= 1; i--){
    vi[k++] = rp[1][i];
  }
  
  return vi;
}
  

// O(ni lg ni + C(R));
vector<int> getpol3(vector<int> vi){
  if(vi.size() <= 6) return getpol(vi);
  set<int> se;
  vector<PT> vw,vc1,vc2,vr;
  vector<int> vi2;
  for(int i = 0; i < vi.size(); i++){
    vw.pb(ve[vi[i]]);
    se.insert(vi[i]);
  }
  vc1 = convex_hull(vw);
  
  for(int i = 0; i < vc1.size(); i++){
    se.erase(vc1[i].i);
  }
  if(se.size() <= 4){
    return getpol(vi);
  }
  
  for(set<int>::iterator it = se.begin(); it != se.end(); it++){
    vr.pb(ve[*it]);
  }
  vc2 = convex_hull(vr);
  
  vector<int> res;
  
  
  
  
  
  
  for(int i = 0; i < vc2.size(); i++){
    vi2.pb(vc2[i].i);
    se.erase(vc2[i].i);
  }
  for(set<int>::iterator it = se.begin(); it != se.end(); it++){
    vi2.pb(ve[*it].i);
  }
  vi2 = getpol1(vi2,0);
  
  
  
  //cout << LinesParallel(ve[vi2[0]], ve[vi2[1]], ve[vi2[1]], ve[vi2[2]]) << endl;
  
  if(!SegmentsIntersect(ve[vi2[0]], vc1.back(), ve[vi2[1]], vc1[0]) || LinesParallel(ve[vi2[0]], ve[vi2[1]], ve[vi2[1]], ve[vi2[2]]) ){
    return getpol(vi);
  }
  
  for(int i = 1; i < vi2.size(); i++){
    res.pb(vi2[i]);
  }
  res.pb(vi2[0]);
  for(int i = 0; i < vc1.size(); i++){
    res.pb(vc1[i].i);
  }
  assert(res.size() == vi.size());
  return res;
}

// O(ni lg ni)
vector<int> getpol(vector<int> vi){
  return getpol1(vi, rand()%vi.size());
}



// O(np lg np)
solution getsol1(){
  vector<int> vi;
  for(int i = 0; i < np; i++){
    vi.pb(i);
  }
  vector<poly> res;
  res.pb(getpol(vi));
  return solution(res);
}



solution localsearch(){
  solution res = getsol1();
  for(int i = 0; i < 2*n; i++){
    res.divide(rand()%res.vp.size());
  }
  return res;
}

solution localsearch2(){
  solution res = getsol1();
  res.divide2();
  return res;
}




class SmallPolygons{
  
  
  public:
  vector <string> choosePolygons(vector <int> points, int N){
    np = points.size();
    n = N;
    for(int i = 0; i < np; i+=2){
      ve.pb(PT(points[i],points[i+1],i/2));
    }
    np/=2;
    
    
    
    int num = 10;
    if(np < 20){
      num = 1000;
    }else if(np < 50){
      num = 500;
    }else if(np < 100){
      num = 200;
    }else if(np < 500){
      num = 100;
    }else{
      num = 50;
    }
    
    solution res = getsol1(), rp;
    /*
    for(int k = 0; k < num; k++){
      rp = localsearch();
      if(rp.getCost() < res.getCost()){
        res = rp;
      }
    }*/
    
    rp = localsearch2();
    if(rp.getCost() < res.getCost()){
      res = rp;
    }
      
    
    
    
    return res.tovs();
  }

};
int main(){
  int Np; 
  vector<int> points;
  int N;
  int v;
  vector<string> ret;
  cin >> Np;
  for (int i=0; i < Np; i++){
    cin >> v;
    points.pb(v);
  }
  cin >> N;
  
  const clock_t begin_time = clock();
  
  
  ret = SmallPolygons().choosePolygons(points, N);
  double dtime = double( clock () - begin_time ) /  CLOCKS_PER_SEC;
  assert(dtime < 30);
  printf("%d\n", ret.size());
  for (int i=0; i < ret.size(); i++){
    printf("%s\n", ret[i].c_str());
    //fprintf(stderr,"VEC: %s\n", ret[i].c_str());
  }
  fflush(stdout);

  return 0;
}
