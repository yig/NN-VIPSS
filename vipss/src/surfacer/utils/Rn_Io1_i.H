// =========================================================
//
//	IO operators on affine1 classes
//
// =========================================================

/** \ingroup SimpleGeometryIO */
//@{
///
template<class Coord>
inline std::ostream& operator<<(std::ostream &out, const R1VectorTC<Coord> &v )
{
    out << v[0] << " ";
    return out;
}

///
template<class Coord>
inline std::istream& operator>>(std::istream &in, R1VectorTC<Coord> &v )
{
    in >> v[0];
    return in;
}

template<class Coord>
inline std::ostream& operator<<(std::ostream &out, const R1CoVectorTC<Coord> &v )
{
    out << v[0] << " ";
    return out;
}

template<class Coord>
inline std::istream& operator>>(std::istream &in, R1CoVectorTC<Coord> &v )
{
    in >> v[0];
    return in;
}


///
template<class Coord>
inline std::ostream& operator<<(std::ostream &out, const R1PointTC<Coord> &p )
{
    out << p[0] << " ";
    return out;
}

///
template<class Coord>
inline std::istream& operator>>(std::istream &in, R1PointTC<Coord> &p )
{
    in >> p[0];
    return in;
}
//@}






