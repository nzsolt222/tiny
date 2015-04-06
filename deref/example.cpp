#include <iostream>
#include <memory>
#include "deref.hpp"

class A
{
public:
    int var;
    void func(){};
};

class B
{
public:
    int var;
    void func(){};
    int& operator*() const noexcept;
};

class C : public B
{
};

struct D
{
    int var;
    void func() const {};
};

struct E
{
    int var;
    void func() const {};
    int& operator*();
};

class F
{
private:
    int var;
    void func() const {};
    int& operator*();
};

template <class T, class U>
void is_same()
{
    static_assert(std::is_same<T, typename nzs::dereference<U>::type>::value,
                  "ERROR");
}

template <class T>
void is_true()
{
    static_assert(nzs::is_dereferancable<T>::value, "ERROR");
}

template <class T>
void is_false()
{
    static_assert(!nzs::is_dereferancable<T>::value, "ERROR");
}

template <class T>
void call_member(T& has_var_and_func)
{
    using namespace nzs;

    deref(has_var_and_func).var = 5;
    deref(has_var_and_func).func();
}

int main()
{
    // is_dereferancable?
    is_true<int**>();
    is_true<std::unique_ptr<int>>();
    is_true<std::unique_ptr<int*>>();
    is_true<std::unique_ptr<int>*>();
    is_true<A*>();
    is_true<B>();
    is_true<C>();
    is_true<D*>();
    is_true<E>();
    is_true<F*>();
    is_true<int (*)()>();

    is_false<int>();
    is_false<A>();
    is_false<D>();
    is_false<F>();

    is_same<int&, int>();
    is_same<int&, int*>();
    is_same<int&, int&>();
    is_same<int*&, int**>();
    is_same<const int&, const int*>();
    is_same<const int*&, const int**>();

    is_same<int&, std::unique_ptr<int>>();
    is_same<int*&, std::unique_ptr<int*>>();
    is_same<std::unique_ptr<int>&, std::unique_ptr<int>*>();

    is_same<A&, A*>();
    is_same<int&, B>();
    is_same<int&, C>();
    is_same<D&, D*>();
    is_same<int&, E>();
    is_same<F&, F*>();

    is_same<int(&)(), int (*)()>();

    A a;
    D d;

    A* pa = new A;
    B* pb = new B;
    C* pc = new C;
    D* pd = new D;
    E* pe = new E;

    std::unique_ptr<A> upa = std::unique_ptr<A>(new A);
    std::unique_ptr<B> upb = std::unique_ptr<B>(new B);
    std::unique_ptr<C> upc = std::unique_ptr<C>(new C);
    std::unique_ptr<D> upd = std::unique_ptr<D>(new D);
    std::unique_ptr<E> upe = std::unique_ptr<E>(new E);

    call_member(a);
    call_member(d);

    call_member(pa);
    call_member(pb);
    call_member(pc);
    call_member(pd);
    call_member(pe);

    call_member(upa);
    call_member(upb);
    call_member(upc);
    call_member(upd);
    call_member(upe);

    delete pa;
    delete pb;
    delete pc;
    delete pd;
    delete pe;

    return 0;
}
