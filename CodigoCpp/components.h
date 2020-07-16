float calculateLocalD(int i,mesh m){
    Matrix matrix;
    Vector row1, row2, row3;

    element e = m.getElement(i);
    
    row1.push_back(calcularTenedor(e, EQUIS, 2, 1, m));
    row1.push_back(calcularTenedor(e, YE, 2, 1, m));
    row1.push_back(calcularTenedor(e, ZETA, 2, 1, m));

    row2.push_back(calcularTenedor(e, EQUIS, 3, 1, m));
    row2.push_back(calcularTenedor(e, YE, 3, 1, m));
    row2.push_back(calcularTenedor(e, ZETA, 3, 1, m));

    row3.push_back(calcularTenedor(e, EQUIS, 4, 1, m));
    row3.push_back(calcularTenedor(e, YE, 4, 1, m));
    row3.push_back(calcularTenedor(e, ZETA, 4, 1, m));

    matrix.push_back(row1);
    matrix.push_back(row2);
    matrix.push_back(row3);

    return determinant(matrix);
}

void calculateAlpha(int i,Matrix &A,mesh m){
    zeroes(A,3);
    element e = m.getElement(i);
    
    A.at(0).at(0) = OperarRestaTenedor(e, YE, ZETA, 3, 4, m);
    A.at(0).at(1) = OperarRestaTenedor(e, YE, ZETA, 4, 2, m);
    A.at(0).at(2) = OperarRestaTenedor(e, YE, ZETA, 2, 3, m);

    A.at(1).at(0) = OperarRestaTenedor(e, EQUIS, ZETA, 4, 3, m);
    A.at(1).at(1) = OperarRestaTenedor(e, EQUIS, ZETA, 2, 4, m);
    A.at(1).at(2) = OperarRestaTenedor(e, EQUIS, ZETA, 3, 2, m);

    A.at(2).at(0) = OperarRestaTenedor(e, EQUIS, YE, 3, 4, m);
    A.at(2).at(1) = OperarRestaTenedor(e, EQUIS, YE, 4, 2, m);
    A.at(2).at(2) = OperarRestaTenedor(e, EQUIS, YE, 2, 3, m);

}

void calculateBeta(Matrix &B){
    zeroes(B,3,4);

    B.at(0).at(0) = -1; 
    B.at(0).at(1) =  1; 
    B.at(0).at(2) =  0; 
    B.at(0).at(3) =  0; 

    B.at(1).at(0) = -1; 
    B.at(1).at(1) = 0; 
    B.at(1).at(2) = 1;
    B.at(1).at(3) = 0;

    B.at(2).at(0) = -1;
    B.at(2).at(1) = 0;
    B.at(2).at(2) = 0;
    B.at(2).at(3) = 1;
}


void calculateGamma(Matrix &m){
	zeroes(m,4,1);
	m.at(0).at(0) = 1;   
	m.at(1).at(0) = 1;  
    m.at(2).at(0) = 1; 
	m.at(3).at(0) = 1;  
}

void calculateAInversa(Matrix &C, mesh m, int i){
    element e = m.getElement(i);
    node n1 = m.getNode(e.getNode1()-1);
    node n2 = m.getNode(e.getNode2()-1);
    node n3 = m.getNode(e.getNode3()-1);
    node n4 = m.getNode(e.getNode4()-1);

    float x1, x2, x3, x4;
    x1 = n1.getX();
    x2 = n2.getX();
    x3 = n3.getX();
    x4 = n4.getX();

    zeroes(C,4,3);
    C.at(0).at(0) = (-(150*pow(x1,3)-59*pow(x1,2)*(x2+x3+x4)+12*x1*pow(x2,2)+x2*(x3+x4)+pow(x3,2)+x3*x4+pow(x4,2))-pow(x2,3)-pow(x2,2)*(x3+x4)-x2*(pow(x3,2)+x3*x4+pow(x4,2))-pow(x3,3)-pow(x3,2)*x4-x3*pow(x4,2)-pow(x4,3))/(840); 
    C.at(1).at(0) = (-(209*pow(x1,3)-71*pow(x1,2)*(2*x2+x3+x4)+13*x1*(3*pow(x2,2)+2*x2*(x3+x4)+pow(x3,2)+x3*x4+pow(x4,2))-4*pow(x2,3)-3*pow(x2,2)*(x3+x4)-2*x2*(pow(x3,2)+x3*x4+pow(x4,2))-pow(x3,3)-pow(x3,2)*x4-x3*pow(x4,2)-pow(x4,3)))/(840); 
    C.at(2).at(0) = (-(209*pow(x1,3)-71*pow(x1,2)*(x2+2*x3+x4)+13*x1*(pow(x2,2)+x2*(2*x3+x4)+3*pow(x3,2)+2*x3*x4+pow(x4,2))-pow(x2,3)-pow(x2,2)*(2*x3+x4)-x2*(3*pow(x3,2)+2*x3*x4+pow(x4,2))-4*pow(x3,3)-3*pow(x3,2)*x4-2*x3*pow(x4,2)-pow(x4,3)))/(840);  
    C.at(3).at(0) = (-(209*pow(x1,3)-71*pow(x1,2)*(x2+x3+2*x4)+13*x1*(pow(x2,2)+x2*(x3+2*x4)+pow(x3,2)+2*x3*x4+3*pow(x4,2))-pow(x2,3)-pow(x2,2)*(x3+2*x4)-x2*(pow(x3,2)+2*x3*x4+3*pow(x4,2))-pow(x3,3)-2*pow(x3,2)*x4-3*x3*pow(x4,2)-4*pow(x4,3)))/(840); 

    C.at(0).at(1) = (-(150*pow(x1,3)-59*pow(x1,2)*(x2+x3+x4)+12*x1*pow(x2,2)+x2*(x3+x4)+pow(x3,2)+x3*x4+pow(x4,2))-pow(x2,3)-pow(x2,2)*(x3+x4)-x2*(pow(x3,2)+x3*x4+pow(x4,2))-pow(x3,3)-pow(x3,2)*x4-x3*pow(x4,2)-pow(x4,3))/(840); 
    C.at(1).at(1) = (-(209*pow(x1,3)-71*pow(x1,2)*(2*x2+x3+x4)+13*x1*(3*pow(x2,2)+2*x2*(x3+x4)+pow(x3,2)+x3*x4+pow(x4,2))-4*pow(x2,3)-3*pow(x2,2)*(x3+x4)-2*x2*(pow(x3,2)+x3*x4+pow(x4,2))-pow(x3,3)-pow(x3,2)*x4-x3*pow(x4,2)-pow(x4,3)))/(840); 
    C.at(2).at(1) = (-(209*pow(x1,3)-71*pow(x1,2)*(x2+2*x3+x4)+13*x1*(pow(x2,2)+x2*(2*x3+x4)+3*pow(x3,2)+2*x3*x4+pow(x4,2))-pow(x2,3)-pow(x2,2)*(2*x3+x4)-x2*(3*pow(x3,2)+2*x3*x4+pow(x4,2))-4*pow(x3,3)-3*pow(x3,2)*x4-2*x3*pow(x4,2)-pow(x4,3)))/(840);  
    C.at(3).at(1) = (-(209*pow(x1,3)-71*pow(x1,2)*(x2+x3+2*x4)+13*x1*(pow(x2,2)+x2*(x3+2*x4)+pow(x3,2)+2*x3*x4+3*pow(x4,2))-pow(x2,3)-pow(x2,2)*(x3+2*x4)-x2*(pow(x3,2)+2*x3*x4+3*pow(x4,2))-pow(x3,3)-2*pow(x3,2)*x4-3*x3*pow(x4,2)-4*pow(x4,3)))/(840); 

    C.at(0).at(2) = (-(150*pow(x1,3)-59*pow(x1,2)*(x2+x3+x4)+12*x1*pow(x2,2)+x2*(x3+x4)+pow(x3,2)+x3*x4+pow(x4,2))-pow(x2,3)-pow(x2,2)*(x3+x4)-x2*(pow(x3,2)+x3*x4+pow(x4,2))-pow(x3,3)-pow(x3,2)*x4-x3*pow(x4,2)-pow(x4,3))/(840); 
    C.at(1).at(2) = (-(209*pow(x1,3)-71*pow(x1,2)*(2*x2+x3+x4)+13*x1*(3*pow(x2,2)+2*x2*(x3+x4)+pow(x3,2)+x3*x4+pow(x4,2))-4*pow(x2,3)-3*pow(x2,2)*(x3+x4)-2*x2*(pow(x3,2)+x3*x4+pow(x4,2))-pow(x3,3)-pow(x3,2)*x4-x3*pow(x4,2)-pow(x4,3)))/(840); 
    C.at(2).at(2) = (-(209*pow(x1,3)-71*pow(x1,2)*(x2+2*x3+x4)+13*x1*(pow(x2,2)+x2*(2*x3+x4)+3*pow(x3,2)+2*x3*x4+pow(x4,2))-pow(x2,3)-pow(x2,2)*(2*x3+x4)-x2*(3*pow(x3,2)+2*x3*x4+pow(x4,2))-4*pow(x3,3)-3*pow(x3,2)*x4-2*x3*pow(x4,2)-pow(x4,3)))/(840);  
    C.at(3).at(2) = (-(209*pow(x1,3)-71*pow(x1,2)*(x2+x3+2*x4)+13*x1*(pow(x2,2)+x2*(x3+2*x4)+pow(x3,2)+2*x3*x4+3*pow(x4,2))-pow(x2,3)-pow(x2,2)*(x3+2*x4)-x2*(pow(x3,2)+2*x3*x4+3*pow(x4,2))-pow(x3,3)-2*pow(x3,2)*x4-3*x3*pow(x4,2)-4*pow(x4,3)))/(840); 

}

float calculatePEscalar(mesh m, int i){
    element e = m.getElement(i);
    node n1 = m.getNode(e.getNode1()-1);
    node n2 = m.getNode(e.getNode2()-1);
    node n3 = m.getNode(e.getNode3()-1);
    node n4 = m.getNode(e.getNode4()-1);

    float y1, y2, y3, y4;
    y1 = n1.getY();
    y2 = n2.getY();
    y3 = n3.getY();
    y4 = n4.getY();
    return ((31*pow(y1,2)-9*y1*(y2+y3+y4)+pow(y2,2)+y2*(y3+y4)+pow(y3,2)+y3*y4+pow(y4,2))/(60));

}

float calculatePe2(mesh m, int i){
     element e = m.getElement(i);
    node n1 = m.getNode(e.getNode1()-1);
    node n2 = m.getNode(e.getNode2()-1);
    node n3 = m.getNode(e.getNode3()-1);
    node n4 = m.getNode(e.getNode4()-1);

    float x1, x2, x3, x4;
    x1 = n1.getX();
    x2 = n2.getX();
    x3 = n3.getX();
    x4 = n4.getX();

    return (-(7*x1-x2-x3-x4))/(24);
}

void calculateOmegaMatrix(Matrix &C, mesh m, int i){
    element e = m.getElement(i);
    node n1 = m.getNode(e.getNode1()-1);
    node n2 = m.getNode(e.getNode2()-1);
    node n3 = m.getNode(e.getNode3()-1);
    node n4 = m.getNode(e.getNode4()-1);

    float y1, y2, y3, y4;
    y1 = n1.getY();
    y2 = n2.getY();
    y3 = n3.getY();
    y4 = n4.getY();

    zeroes(C,4,3);
    C.at(0).at(0) = ((39*pow(y1,2)-10*y1*(y2+y3+y4)+pow(y2,2)+y2*(y3+y4)+pow(y3,2)+y3*y4+pow(y4,2))/(360)); 
    C.at(1).at(0) = ((49*pow(y1,2)-11*y1*(2*y2+y3+y4)+3*pow(y2,2)+2*y2*(y3+y4)+pow(y3,2)+y3*y4+pow(y4,2))/(360)); 
    C.at(2).at(0) = ((49*pow(y1,2)-11*y1*(y2+2*y3+y4)+pow(y2,2)+y2*(2*y3+y4)+3*pow(y3,2)+2*y3*y4+pow(y4,2))/(360));  
    C.at(3).at(0) = ((49*pow(y1,2)-11*y1*(y2+y3+2*y4)+pow(y2,2)+y2*(y3+2*y4)+pow(y3,2)+2*y3*y4+3*pow(y4,2))/(360)); 

    C.at(0).at(1) = ((39*pow(y1,2)-10*y1*(y2+y3+y4)+pow(y2,2)+y2*(y3+y4)+pow(y3,2)+y3*y4+pow(y4,2))/(360)); 
    C.at(1).at(1) = ((49*pow(y1,2)-11*y1*(2*y2+y3+y4)+3*pow(y2,2)+2*y2*(y3+y4)+pow(y3,2)+y3*y4+pow(y4,2))/(360)); 
    C.at(2).at(1) = ((49*pow(y1,2)-11*y1*(y2+2*y3+y4)+pow(y2,2)+y2*(2*y3+y4)+3*pow(y3,2)+2*y3*y4+pow(y4,2))/(360));  
    C.at(3).at(1) = ((49*pow(y1,2)-11*y1*(y2+y3+2*y4)+pow(y2,2)+y2*(y3+2*y4)+pow(y3,2)+2*y3*y4+3*pow(y4,2))/(360)); 

    C.at(0).at(2) = ((39*pow(y1,2)-10*y1*(y2+y3+y4)+pow(y2,2)+y2*(y3+y4)+pow(y3,2)+y3*y4+pow(y4,2))/(360)); 
    C.at(1).at(2) = ((49*pow(y1,2)-11*y1*(2*y2+y3+y4)+3*pow(y2,2)+2*y2*(y3+y4)+pow(y3,2)+y3*y4+pow(y4,2))/(360)); 
    C.at(2).at(2) = ((49*pow(y1,2)-11*y1*(y2+2*y3+y4)+pow(y2,2)+y2*(2*y3+y4)+3*pow(y3,2)+2*y3*y4+pow(y4,2))/(360));  
    C.at(3).at(2) = ((49*pow(y1,2)-11*y1*(y2+y3+2*y4)+pow(y2,2)+y2*(y3+2*y4)+pow(y3,2)+2*y3*y4+3*pow(y4,2))/(360)); 

}



float calculateLocalJ(int i,mesh m){
    Matrix matrix;
    Vector row1, row2, row3;

    element e = m.getElement(i);
    
    row1.push_back(calcularTenedor(e, EQUIS, 2, 1, m));
    row1.push_back(calcularTenedor(e, EQUIS, 3, 1, m));
    row1.push_back(calcularTenedor(e, EQUIS, 4, 1, m));

    row2.push_back(calcularTenedor(e, YE, 2, 1, m));
    row2.push_back(calcularTenedor(e, YE, 3, 1, m));
    row2.push_back(calcularTenedor(e, YE, 4, 1, m));

    row3.push_back(calcularTenedor(e, ZETA, 2, 1, m));
    row3.push_back(calcularTenedor(e, ZETA, 3, 1, m));
    row3.push_back(calcularTenedor(e, ZETA, 4, 1, m));

    matrix.push_back(row1);
    matrix.push_back(row2);
    matrix.push_back(row3);

    return determinant(matrix);
}

Matrix createLocalM(int e,mesh &m){
    Matrix matrixAInversa,matrixE,matrixK,matrixL,matrixf,matrixh,matrixC;
    float const1,const2,J,Determinant;

    /* [ C  -K ]
       [ E  -L ]
    */

    //Matrix C
    Matrix Alpha, Beta;
    Determinant = calculateLocalD(e,m);
    J = calculateLocalJ(e,m);
    
    if(Determinant == 0){
        cout << "\n!---CATASTROPHIC FAILURE---M!\n";
        exit(EXIT_FAILURE);
    }
    float real_a = (float) (J)/(Determinant);
    calculateAInversa(matrixAInversa,m,e);
    calculateAlpha(e,Alpha,m);
    calculateBeta(Beta);
    productRealMatrix(real_a, productMatrixMatrix(matrixAInversa,productMatrixMatrix(Alpha,Beta,3,3,4),4,3,4),matrixC);

    //Matrix K
    Matrix Alpha_t,Beta_t;
    float p = calculatePEscalar(m,e);    
    float real_k = (float) (p*J)/(Determinant*Determinant);
    transpose(Alpha,Alpha_t);
    transpose(Beta,Beta_t);
    productRealMatrix(-real_k,productMatrixMatrix(Beta_t,productMatrixMatrix(Alpha_t,productMatrixMatrix(Alpha,Beta,3,3,4),3,3,4),4,3,4),matrixK);

    //Matrix E
    Matrix Omega;
    float real_g = (float) (J)/(Determinant);
    calculateOmegaMatrix(Omega,m,e);
    productRealMatrix(real_g,productMatrixMatrix(Omega,productMatrixMatrix(Alpha,Beta,3,3,4),4,3,4),matrixE);

    //Matrix L
    Matrix Omega_t;
    float p2 = calculatePe2(m,e);
    float real_d = (float)(J/(720*Determinant));
    productRealMatrix(-real_d,productMatrixMatrix(Beta_t,productMatrixMatrix(Alpha_t,productMatrixMatrix(Alpha,Beta,3,3,4),3,3,4),4,3,4),matrixL);

    //Matrix M
    Matrix M;
    zeroes(M,8);
    ubicarSubMatriz(M,0,3,0,3, matrixC);
    ubicarSubMatriz(M,0,3,4,7,matrixK);
    ubicarSubMatriz(M,4,7,0,3,matrixE);
    ubicarSubMatriz(M,4,7,4,7,matrixL);
    return M;
}

Vector createLocalb(int e,mesh &m){

    float J;
    Vector b;
    J = calculateLocalJ(e,m);

    if(J == 0){
        cout << "\n!---CATASTROPHIC FAILURE---J!\n";
        exit(EXIT_FAILURE);
    }
    zeroes(b,8);
    b.at(0) = (J*28.4)/24;
    b.at(1) = (J*28.4)/24;
    b.at(2) = (J*28.4)/24;
    b.at(3) = (J*28.4)/24;

    b.at(4) = (J*52.8)/24;
    b.at(5) = (J*52.8)/24;
    b.at(6) = (J*52.8)/24;
    b.at(7) = (J*52.8)/24;

    return b;
}
