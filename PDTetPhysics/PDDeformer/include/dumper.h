#pragma once

#ifndef NO_PHYSBAM
// #include "Discretization.h"
#include "Geometry.h"
#include "Iterator.h"
#endif

#include <iostream>
#include <fstream>
#include <string>
#include <vector>

namespace dumper {

    template <class T, class IntType>
        void writeDense(const IntType n,  const T* a, const std::string fileName = "dense_matrix.txt") {
        std::ofstream ofs;
        ofs.open(fileName);
        ofs<<"a = [";

        // print in column major
        for (int j=0; j<n; j++) {
            for (int i=0; i<n; i++)
                ofs<<a[j*n+i]<<" ";
            ofs<<";";
        }
        ofs<<"]";
        ofs.close();
    }

    // out put to matlab sparse format
    template <class T, class IntType>
        void writeSparse(const IntType n, const IntType* ia, const IntType* ja, const T* a, const std::string fileName = "matrix.txt") {
        std::ofstream ofs;
        ofs.open(fileName);
        ofs<<"i = [";
        for (int i=0; i<n; i++)
            for (int jj=ia[i]; jj<ia[i+1]; jj++)
                ofs<<i<<" ";
        ofs<<"];"<<std::endl;

        ofs<<"j = [";
        for (int i=0; i<n; i++)
            for (int jj=ia[i]; jj<ia[i+1]; jj++)
                ofs<<ja[jj]<<" ";
        ofs<<"];"<<std::endl;

        ofs<<"a = [";
        for (int i=0; i<n; i++)
            for (int jj=ia[i]; jj<ia[i+1]; jj++)
                ofs<<a[jj]<<" ";
        ofs<<"];"<<std::endl;

        ofs.close();
    }

    template <class T, class IntType>
        void writeByte(const IntType n, const T* a, std::ofstream& ofs) {
        if ( ofs )
            ofs.write( reinterpret_cast<const char*>( a ), sizeof(T)*n);
    }

    template <class T, class IntType>
        void writeCSRbyte(const IntType n, const IntType* ia, const IntType* ja, const T* a, const IntType m=0, const std::string iFileName = "i.txt", const std::string aFileName = "a.txt") {
// dump data to files
        std::ofstream ofs;
        ofs.open( iFileName, std::ios::binary );
        if ( ofs )  {
            ofs.write( reinterpret_cast<const char*>( &n ), sizeof(IntType));
            ofs.write( reinterpret_cast<const char*>( ia ), sizeof(IntType)*(n+1) );
            ofs.write( reinterpret_cast<const char*>( ja ), sizeof(IntType)*(ia[n]));
            ofs.write( reinterpret_cast<const char*>( &m ), sizeof(IntType));
            // Close the file to unlock it
            ofs.close();
            ofs.open( aFileName, std::ios::binary );
            if ( ofs )  {
                ofs.write( reinterpret_cast<const char*>( a ), sizeof(T)*(ia[n]));
                ofs.close();
            }
        }

    }

    template <class T, class IntType>
        void readByte(const IntType n, T* a, std::ifstream& ifs) {
        if ( ifs )
            ifs.read( reinterpret_cast<char*>( a ), sizeof(T)*n);
    }


    template <class T, class IntType>
        void readCSRbyte(IntType* n,  IntType** ia,  IntType** ja,  T** a, IntType* m, const std::string iFileName = "i.txt", const std::string aFileName = "a.txt") {
        // Use a new object so we don't have to worry
        // about error states in the old object
        std::ifstream ifs;
        ifs.open( iFileName, std::ios::binary );
        if ( ifs )  {
            ifs.read( reinterpret_cast<char*>( n ), sizeof(IntType));
            std::cout<<"n = "<<*n<<std::endl;
            *ia = new IntType[*n+1];
            ifs.read( reinterpret_cast<char*>( *ia ), sizeof(IntType)*(*n+1) );
            std::cout<<"read in ia"<<std::endl;
            IntType nnz = (*ia)[*n];
            std::cout<<"nnz = "<<nnz<<std::endl;
            *ja = new IntType[nnz];
            ifs.read( reinterpret_cast<char*>( *ja ), sizeof(IntType)*nnz);
            std::cout<<"read in ja"<<std::endl;
            ifs.read( reinterpret_cast<char*>( m ), sizeof(IntType));
            std::cout<<"m = "<<*m<<std::endl;
            // Close the file to unlock it
            ifs.close();
            ifs.open( aFileName, std::ios::binary );
            if ( ifs )  {
                *a = new T[nnz];
                ifs.read( reinterpret_cast<char*>( *a ), sizeof(T)*nnz);
                std::cout<<"read in a"<<std::endl;
                ifs.close();
            }


        }
    }

#ifndef NO_PHYSBAM
    template<int ElementNodes> void writeElementsByte(const std::vector<std::array<int, ElementNodes>> &elements) {
        std::cout << "elements.size() = " << elements.size() << std::endl;
        std::ofstream fs;
        fs.open("elements.dat", std::ios::out | std::ios::binary);
        int nElements = elements.size();
        fs.write( reinterpret_cast<const char*>(&nElements), sizeof(int) );
        fs.write( reinterpret_cast<const char*>(&elements[0][0]), sizeof(int) * nElements * ElementNodes );
        fs.close();
    }

    template<class T, int d> void writePositionsByte(const std::vector<PhysBAM::VECTOR<T, d>> &positions) {
        std::cout << "positions.size() = " << positions.size() << std::endl;
        std::ofstream fs;
        fs.open("positions.dat", std::ios::out | std::ios::binary);
        int nParticles = positions.size();
        fs.write( reinterpret_cast<const char*>(&nParticles), sizeof(int) );
        for(const auto& x : positions)
            for(int v = 1; v <= d; v++){
                float xFloat = float(x(v));
                fs.write( reinterpret_cast<const char*>( &xFloat ), sizeof(float) );
            }
        fs.close();
    }

    template<int ElementNodes> void writeElements(const std::vector<std::array<int, ElementNodes>> &elements) {
        std::ofstream Myout("Elements.txt");
        for (int i = 0; i < elements.size(); i++) {
            const auto &element = elements[i];
            for (int j = 0; j < ElementNodes; j++)
                Myout << element[j] << " ";
            Myout << std::endl;
        }
    }

    template<int ElementNodes> void readElements(std::vector<std::array<int, ElementNodes>> &elements, const std::string prefix="") {
        std::ifstream Myin(prefix+"Elements.txt");
        std::array<int, ElementNodes> element;
        std::string line;
        while (std::getline(Myin, line)) {
            std::istringstream iss(line);
            int index;
            for (int i = 0; i < ElementNodes; i++) {
                iss >> index;
                element[i] = index;
            }
            elements.push_back(element);
        }
    }


    template<class T, int d> void writePositions(const std::vector<PhysBAM::VECTOR<T, d>> &positions) {
        std::ofstream Myout("Positions.txt");
        for (int i = 0; i < positions.size(); i++) {
            const auto &position = positions[i];
            for (int j = 0; j < d; j++)
                Myout << position(j + 1) << " ";
            Myout << std::endl;
        }
    }

    template<class T, int d> void readPositions(std::vector<PhysBAM::VECTOR<T, d>> &positions, const std::string prefix="") {
        std::ifstream Myin(prefix+"Positions.txt");
        PhysBAM::VECTOR<T, d> position;
        std::string line;
        while (std::getline(Myin, line)) {
            std::istringstream iss(line);
            T x;
            for (int i = 0; i < d; i++) {
                iss >> x;
                position(i + 1) = x;
            }
            positions.push_back(position);
        }
    }

    template<int d> void writeNodeTypes(const std::vector<typename PhysBAM::Geometry<d>::NodeType>  &nodeTypes) {
        std::ofstream Myout("NodeTypes.txt");
        for (int i = 0; i < nodeTypes.size(); i++)
            Myout << static_cast<int>(nodeTypes[i]) << std::endl;
    }

    template<int d> void readNodeTypes(std::vector<typename PhysBAM::Geometry<d>::NodeType>  &nodeTypes, const std::string prefix="") {
        std::ifstream Myin(prefix+"NodeTypes.txt");
        std::string line;
        while (std::getline(Myin, line)) {
            std::istringstream iss(line);
            int nodeType;
            iss >> nodeType;
            nodeTypes.push_back(static_cast<typename PhysBAM::Geometry<d>::NodeType>(nodeType));
        }
    }

    template<class T, int d, int ElementNodes> void writeConstraints(const std::vector<PhysBAM::SoftConstraint<PhysBAM::VECTOR<T, d>, ElementNodes, int>> &constraints) {
        std::ofstream Myout("Constraints.txt");
        for (int i = 0; i < constraints.size(); i++) {
            const auto &constraint = constraints[i];
            for (int j = 0; j < ElementNodes; j++)
                Myout << constraint.m_elementIndex[j] << " ";
            for (int j = 0; j < d; j++)
                Myout << constraint.m_xT(j + 1) << " ";
            Myout << constraint.m_stiffness << " ";
            for (int j = 0; j < d + 1; j++)
                Myout << constraint.m_weights[j] << " ";
            Myout << std::endl;
        }
    }

    template<class T, int d, int ElementNodes> void writeSutures(const std::vector<PhysBAM::SutureConstraint<PhysBAM::VECTOR<T, d>, ElementNodes, int>>& constraints) {
        std::ofstream Myout("Sutures.txt");
        for (int i = 0; i < constraints.size(); i++) {
            const auto& constraint = constraints[i];
            for (int j = 0; j < ElementNodes; j++)
                Myout << constraint.m_elementIndex1[j] << " ";
            for (int j = 0; j < ElementNodes; j++)
                Myout << constraint.m_elementIndex2[j] << " ";
            Myout << constraint.m_stiffness << " ";
            for (int j = 0; j < d + 1; j++)
                Myout << constraint.m_weights1[j] << " ";
            for (int j = 0; j < d + 1; j++)
                Myout << constraint.m_weights2[j] << " ";
            Myout << std::endl;
        }
    }

    template<class T, int d, int ElementNodes> void readConstraints(std::vector<PhysBAM::SoftConstraint<PhysBAM::VECTOR<T, d>, ElementNodes, int>> &constraints, const std::string prefix="") {
        std::ifstream Myin(prefix+"Constraints.txt");
        PhysBAM::SoftConstraint<PhysBAM::VECTOR<T, d>, ElementNodes, int> constraint;
        std::string line;
        while (std::getline(Myin, line)) {
            std::istringstream iss(line);
            for (int i = 0; i < ElementNodes; i++)
                iss >> constraint.m_elementIndex[i];
            for (int i = 0; i < d; i++)
                iss >> constraint.m_xT(i + 1);
            iss >> constraint.m_stiffness;
            for (int i = 0; i < d + 1; i++)
                iss >> constraint.m_weights[i];
            constraints.push_back(constraint);
        }
    }

    template<class T, int d, int ElementNodes> void writeCollisionConstraints(const std::vector<PhysBAM::SoftConstraint<PhysBAM::VECTOR<T, d>, ElementNodes, int>> &constraints) {
        std::ofstream Myout("CollisionConstraints.txt");
        for (int i = 0; i < constraints.size(); i++) {
            const auto &constraint = constraints[i];
            for (int j = 0; j < ElementNodes; j++)
                Myout << constraint.m_elementIndex[j] << " ";
            for (int j = 0; j < d; j++)
                Myout << constraint.m_xT(j + 1) << " ";
            Myout << constraint.m_stiffness << " ";
            for (int j = 0; j < d + 1; j++)
                Myout << constraint.m_weights[j] << " ";
            Myout << std::endl;
        }
    }

    template<class T, int d, int ElementNodes> void readCollisionConstraints(std::vector<PhysBAM::SoftConstraint<PhysBAM::VECTOR<T, d>, ElementNodes, int>> &constraints, const std::string prefix="") {
        std::ifstream Myin(prefix+"CollisionConstraints.txt");
        PhysBAM::SoftConstraint<PhysBAM::VECTOR<T, d>, ElementNodes, int> constraint;
        std::string line;
        while (std::getline(Myin, line)) {
            std::istringstream iss(line);
            for (int i = 0; i < ElementNodes; i++)
                iss >> constraint.m_elementIndex[i];
            for (int i = 0; i < d; i++)
                iss >> constraint.m_xT(i + 1);
            iss >> constraint.m_stiffness;
            for (int i = 0; i < d + 1; i++)
                iss >> constraint.m_weights[i];
            constraints.push_back(constraint);
        }
    }
    #endif
}
