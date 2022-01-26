#pragma once

#include "Iterator.h"

namespace PhysBAM {

    template <class StateVariableType> struct Algebra;

    template <class StateVariable> struct Algebra {
        using StateVariableType = StateVariable;
        using IteratorType = Iterator<StateVariableType>;
        using VectorType = typename IteratorType::DataType;
        using T = typename VectorType::ELEMENT;

        static T innerProduct(const StateVariableType &v, const StateVariableType &w) {
            T result = 0.;
            for (IteratorType iterator(v); !iterator.isEnd(); iterator.next())
                result += VectorType::Dot_Product(iterator.value(v), iterator.value(w));
            return result;
        }

        static T infinityNorm(const StateVariableType &v) {
            // this is a weird kind of infinity norm...
            T result = 0.;
            for (IteratorType iterator(v); !iterator.isEnd(); iterator.next())
                result = std::max(result, iterator.value(v).Magnitude());
            return result;
        }

        static T lNorm(const StateVariableType &v, T l) {
            // this is a weird kind of l norm...
            T result = 0.;
            for (IteratorType iterator(v); !iterator.isEnd(); iterator.next())
                result += pow(iterator.value(v).Magnitude(), l);
            return pow(result, 1/l);
        }

        static void negate(StateVariableType &v) {
            for (IteratorType iterator(v); !iterator.isEnd(); iterator.next())
                iterator.value(v) *= T(-1);
        }

        static void clear(StateVariableType &v) {
            IteratorType iterator(v);
            iterator.resize(v);
        }

        static void addTo(StateVariableType &X, const StateVariableType &dX) {
            for (IteratorType iterator(X); !iterator.isEnd(); iterator.next())
                iterator.value(X) += iterator.value(dX);
        }

        static void subtractFrom(StateVariableType &X, const StateVariableType &dX) {
            for (IteratorType iterator(X); !iterator.isEnd(); iterator.next())
                iterator.value(X) -= iterator.value(dX);
        }

        static void multiplyBy(StateVariableType &X, const T alpha) {
            for (IteratorType iterator(X); !iterator.isEnd(); iterator.next())
                iterator.value(X) *= alpha;
        }

        static void scaleAndCopy(StateVariableType &Y, const StateVariableType &X, const T alpha) {
            for (IteratorType iterator(X); !iterator.isEnd(); iterator.next())
                iterator.value(Y) = alpha * iterator.value(X);
        }

        static void axpyAndCopy(StateVariableType &Z, const StateVariableType &X, const StateVariableType &B,
                                const T alpha) {
            for (IteratorType iterator(X); !iterator.isEnd(); iterator.next())
                iterator.value(Z) = alpha * iterator.value(X) + iterator.value(B);
        }



    };
}
