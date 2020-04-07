/*__global__ void getPointEvals(float* unknowns, float* mPoints, float* outs)
{
    int mpidx = i * 3;
    int cidx = i * 4;
    Vector3f vk((*mPoints)[mpidx], (*mPoints)[mpidx + 1], (*mPoints)[mpidx + 2]);
    float alpha = (*unknowns)(cidx);
    Vector3f beta((*unknowns)(cidx + 1), (*unknowns)(cidx + 2), (*unknowns)(cidx + 3));
    Vector3f diff = p - vk;
    Vector3f grad(derivx(diff(0), diff(1), diff(2)), derivy(diff(0), diff(1), diff(2)), derivz(diff(0), diff(1), diff(2)));
    //std::cout << alpha << std::endl;
    //std::cout << beta(0) << "," << beta(1) << "," << beta(2) << " " << diff(0) << "," << diff(1) << "," << diff(2) << " " << grad(0) << "," << grad(1) << "," << grad(2) << std::endl;
    out += alpha * smoothfunc(diff(0), diff(1), diff(2)) - beta.dot(grad);
    //std::cout << out << std::endl;
}*/