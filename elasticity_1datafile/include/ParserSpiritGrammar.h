//@HEADER
/*
*******************************************************************************

    Copyright (C) 2004, 2005, 2007 EPFL, Politecnico di Milano, INRIA
    Copyright (C) 2010 EPFL, Politecnico di Milano, Emory University

    This file is part of LifeV.

    LifeV is free software; you can redistribute it and/or modify
    it under the terms of the GNU Lesser General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    LifeV is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    Lesser General Public License for more details.

    You should have received a copy of the GNU Lesser General Public License
    along with LifeV.  If not, see <http://www.gnu.org/licenses/>.

*******************************************************************************
*/
//@HEADER

/*!
 *  @file
 *  @brief File containing the Parser grammar
 *
 *  @date 05-02-2010
 *  @author Cristiano Malossi <cristiano.malossi@epfl.ch>
 *
 *  @maintainer Cristiano Malossi <cristiano.malossi@epfl.ch>
 */

#ifndef Parser_SpiritGrammar_H
#define Parser_SpiritGrammar_H 

#include "ParserDefinitions.h"

namespace LifeV
{

//#ifdef HAVE_BOOST_SPIRIT_QI

//! ParserSpiritGrammar - A string parser grammar based on \c boost::spirit::qi
/*!
 *  @author(s) Cristiano Malossi
 *
 *  \c ParserSpiritGrammar is a \c boost::spirit::qi based class to perform
 *  evaluation of \c std::string expressions.
 *
 *  <b>EXAMPLE - HOW TO USE</b>
 *  Let's consider the following example: suppose that we have this function:
 *
 *  [u,v,w] = f(x,y,z,t)
 *
 *  where
 *
 *  u(x)   = a*b*x
 *  v(x,y) = a/b*sqrt(x^2 + y^2)
 *  w(t)   = b*t;
 *
 *  with "a" and "b" constants such that a=5.12345, b=9.999999.
 *
 *  To evaluate function f(x,y,z,t), we use this syntax:
 *
 *  string = "a=5.12345 ; b=9.999999 ; (a*b*x, a/b*sqrt(x^2 + y^2), b*t)"
 *
 *  where semicolons (";") separate constants and commas (",") separate output functions.
 *
 *  NOTE:
 *  Currently ParserSpiritGrammar works with the following operators:
 *  \verbatim
 *  +, -, *, /, ^, sqrt(), sin(), cos(), tan(), exp(), log(), log10(), >, <.
 *  \endverbatim
 *
 */
template < typename IteratorType = std::string::const_iterator, typename ResultsType = std::vector < scalar_type > >
class ParserSpiritGrammar : public qi::grammar< IteratorType, ResultsType(), ascii::space_type >
{
public:

    //! @name Public Types
    //@{

    typedef IteratorType                                        iterator_Type;
    typedef boost::iterator_range< iterator_Type >              iteratorRange_Type;
    typedef ResultsType                                         results_Type;

    //@}


    //! @name Constructors & Destructor
    //@{

    //! Constructor
    explicit ParserSpiritGrammar();

    //! Copy constructor
    /*!
     * @param ParserSpiritGrammar ParserSpiritGrammar
     */
    explicit ParserSpiritGrammar( const ParserSpiritGrammar& spiritGrammar );

    //! Destructor
    virtual ~ParserSpiritGrammar() {}

    //@}


    //! @name Operators
    //@{

    //! Operator =
    /*!
     * @param SpiritGrammar ParserSpiritGrammar
     * @return reference to a copy of the class
     */
    ParserSpiritGrammar& operator=( const ParserSpiritGrammar& spiritGrammar );

    //@}


    //! @name Methods
    //@{

    /*! Assign a variable using a \c boost::iterator_range
     *
     * @param stringIteratorRange name of the parameter
     * @param value value of the parameter
     */
    void assignVariable( const iteratorRange_Type& stringIteratorRange, const scalar_type& value ) { setVariable( std::string( stringIteratorRange.begin(), stringIteratorRange.end() ), value ); }

    //! Clear all the variables.
    void clearVariables() { M_variable.clear(); }

    /*
        //! Show all the variables
        void ShowMe() { std::cout << "ShowMe called!" << std::endl; }
    */

    //@}


    //! @name Set Methods
    //@{

    //! Set default variables
    void setDefaultVariables();

    /*! Set/replace a variable
     *
     * @param name name of the parameter
     * @param value value of the parameter
     */
    void setVariable( const std::string& name, const scalar_type& value );

    //@}

    //! @name Get Methods
    //@{

    /*! Get variable
     *
     * @param name name of the parameter
     * @return value of the variable
     */
    scalar_type& variable( const std::string& name ) { return M_variable.at( name ); }

    //@}

private:

    //! @name Private Methods
    //@

    /*! Phoenix wrapper for \c std::sin
     *
     * @param value input value
     * @return sin( value )
     */
    scalar_type sin( const scalar_type& value ) const { return std::sin( value ); }

    /*! Phoenix wrapper for \c std::cos
     *
     * @param value input value
     * @return cos( value )
     */
    scalar_type cos( const scalar_type& value ) const { return std::cos( value ); }

    /*! Phoenix wrapper for \c std::tan
     *
     * @param value input value
     * @return tan( value )
     */
    scalar_type tan( const scalar_type& value ) const { return std::tan( value ); }

    /*! Phoenix wrapper for \c std::pow
     *
     * @param base input base
     * @param exponent input exponent
     * @return pow( base, exponent )
     */
    scalar_type pow( const scalar_type& Base, const scalar_type& Exponent ) const { return std::pow( Base, Exponent ); }

    /*! Phoenix wrapper for \c std::sqrt
     *
     * @param value input value
     * @return sqrt( value )
     */
    scalar_type sqrt( const scalar_type& value ) const { return std::sqrt( value ); }

    /*! Phoenix wrapper for \c std::exp
     *
     * @param value input value
     * @return exp( value )
     */
    scalar_type exp( const scalar_type& value ) const { return std::exp( value ); }

    /*! Phoenix wrapper for \c std::log
     *
     * @param value input value
     * @return log( value )
     */
    scalar_type log( const scalar_type& value ) const { return std::log( value ); }

    /*! Phoenix wrapper for \c std::log10
     *
     * @param value input value
     * @return log10( value )
     */
    scalar_type log10( const scalar_type& value ) const { return std::log10( value ); }

    //@}

    qi::rule< iterator_Type, results_Type(), ascii::space_type > M_start;

    qi::rule< iterator_Type, void(), ascii::space_type >         M_assignment;
//    qi::rule< iterator_Type, void(), ascii::space_type >         M_command;

    qi::rule< iterator_Type, double(), ascii::space_type >       M_expression;
    qi::rule< iterator_Type, double(), ascii::space_type >       M_compare;
    qi::rule< iterator_Type, double(), ascii::space_type >       M_plusMinus;
    qi::rule< iterator_Type, double(), ascii::space_type >       M_multiplyDivide;
    qi::rule< iterator_Type, double(), ascii::space_type >       M_elevate;
    qi::rule< iterator_Type, double(), ascii::space_type >       M_element;
    qi::rule< iterator_Type, double(), ascii::space_type >       M_number;
    qi::rule< iterator_Type, double(), ascii::space_type >       M_function;
    qi::rule< iterator_Type, double(), ascii::space_type >       M_group;

    qi::symbols<char, scalar_type >                                     M_variable;
};



// ===================================================
// Constructors & Destructor
// ===================================================
template < typename IteratorType, typename ResultsType >
ParserSpiritGrammar< IteratorType, ResultsType >::ParserSpiritGrammar() :
        ParserSpiritGrammar::base_type( M_start ),
        M_start                               (),
        M_assignment                          (),
//        M_command                             (),
        M_expression                          (),
        M_compare                             (),
        M_plusMinus                           (),
        M_multiplyDivide                      (),
        M_elevate                             (),
        M_element                             (),
        M_number                              (),
        M_function                            (),
        M_group                               (),
        M_variable                            ()
{
    M_start
    =
        (
            M_assignment
//        |  M_command
        |  ( -qi::lit('[') >> M_expression % ',' >> -qi::lit(']') )
        )
        ;

    M_assignment
    =
        (    qi::raw[qi::lexeme[(qi::alpha | '_') >> *(qi::alnum | '_')]]
        >>   qi::lit('=')
        >>   M_expression
        )              [phoenix::bind(&ParserSpiritGrammar::assignVariable,this, qi::_1, qi::_2)]
        ;
    /*
        M_command
            =  qi::lit("ShowMe")[phoenix::bind(&ParserSpiritGrammar::ShowMe, this)]
        ;
    */

    M_expression
    =  *M_compare                                        [qi::_val = qi::_1]
       ;

    M_compare
    =   M_plusMinus                                      [qi::_val = qi::_1]
        >> *(
  	    qi::lit('>') >> M_plusMinus                  [qi::_val = qi::_val > qi::_1]
        |   qi::lit('?') >> M_plusMinus                  [qi::_val = (qi::_val == qi::_1)]
        |   qi::lit('<') >> M_plusMinus                  [qi::_val = qi::_val < qi::_1]
        )
        ;

    M_plusMinus
    =   M_multiplyDivide                                 [qi::_val = qi::_1]
        >> *(
            qi::lit('+') >> M_multiplyDivide             [qi::_val += qi::_1]
        |   qi::lit('-') >> M_multiplyDivide             [qi::_val -= qi::_1]
        )
        ;

    M_multiplyDivide
    =   M_elevate                                        [qi::_val = qi::_1]
        >> *(
            qi::lit('*') >> M_elevate                    [qi::_val *= qi::_1]
        |   qi::lit('/') >> M_elevate                    [qi::_val /= qi::_1]
        )
        ;

    M_elevate
    =
        (
            qi::lit('-') >> M_element                    [qi::_val = qi::_1]
            >>  (
                qi::lit('^') >> M_element                [qi::_val = -phoenix::bind(&ParserSpiritGrammar::pow, this, qi::_val, qi::_1)]
            )
            >> *(
                qi::lit('^') >> M_element                [qi::_val = phoenix::bind(&ParserSpiritGrammar::pow, this, qi::_val, qi::_1)]
            )
        )
        |
        (
            M_element                                    [qi::_val = qi::_1]
            >> *(
                qi::lit('^') >> M_element                [qi::_val = phoenix::bind(&ParserSpiritGrammar::pow, this, qi::_val, qi::_1)]
            )
        )
        ;

    M_element
    =
        (
            qi::lit('-') >> M_element                    [qi::_val = -qi::_1]
        )
        |
        (
            M_number                                     [qi::_val = qi::_1]
        |   M_function                                   [qi::_val = qi::_1]
        |   M_variable                                   [qi::_val = qi::_1]
        |   M_group                                      [qi::_val = qi::_1]
        )
        ;

    M_number
    =
        (
            qi::double_
//            ||   ('.' >> qi::double_)
//            >>   -('.' >> qi::double_) | ('.' >> qi::double_)
        )
        ;

    M_function
    =
        (
            qi::lit("sin")   >> M_group                  [qi::_val = phoenix::bind(&ParserSpiritGrammar::sin, this,   qi::_1)]
        |   qi::lit("cos")   >> M_group                  [qi::_val = phoenix::bind(&ParserSpiritGrammar::cos, this,   qi::_1)]
        |   qi::lit("tan")   >> M_group                  [qi::_val = phoenix::bind(&ParserSpiritGrammar::tan, this,   qi::_1)]
        |   qi::lit("sqrt")  >> M_group                  [qi::_val = phoenix::bind(&ParserSpiritGrammar::sqrt, this,  qi::_1)]
        |   qi::lit("exp")   >> M_group                  [qi::_val = phoenix::bind(&ParserSpiritGrammar::exp, this,   qi::_1)]
        |   qi::lit("log")   >> M_group                  [qi::_val = phoenix::bind(&ParserSpiritGrammar::log, this,   qi::_1)]
        |   qi::lit("log10") >> M_group                  [qi::_val = phoenix::bind(&ParserSpiritGrammar::log10, this, qi::_1)]
        )
        ;

    M_group
    =
        (
            '('
        >>  M_expression                                 [qi::_val = qi::_1]
        >>  ')'
        )
        ;
}

template < typename IteratorType, typename ResultsType >
ParserSpiritGrammar< IteratorType, ResultsType >::ParserSpiritGrammar( const ParserSpiritGrammar& spiritGrammar ) :
        ParserSpiritGrammar::base_type     ( spiritGrammar.M_start ),
        M_start                               ( spiritGrammar.M_start ),
        M_assignment                          ( spiritGrammar.M_assignment ),
//        M_command                             ( spiritGrammar.M_command ),
        M_expression                          ( spiritGrammar.M_expression ),
        M_compare                             ( spiritGrammar.M_compare ),
        M_plusMinus                           ( spiritGrammar.M_plusMinus ),
        M_multiplyDivide                      ( spiritGrammar.M_multiplyDivide ),
        M_elevate                             ( spiritGrammar.M_elevate ),
        M_element                             ( spiritGrammar.M_element ),
        M_number                              ( spiritGrammar.M_number ),
        M_function                            ( spiritGrammar.M_function ),
        M_group                               ( spiritGrammar.M_group ),
        M_variable                            ( spiritGrammar.M_variable )
{
}

// ===================================================
// Operators
// ===================================================
template < typename IteratorType, typename ResultsType >
ParserSpiritGrammar< IteratorType, ResultsType >&
ParserSpiritGrammar< IteratorType, ResultsType >::operator=( const ParserSpiritGrammar& spiritGrammar )
{
    if ( this != &spiritGrammar )
    {
        ParserSpiritGrammar::base_type::operator=( spiritGrammar.M_start );
//        ParserSpiritGrammar< IteratorType, ResultsType >::operator=( spiritGrammar );
        M_start                               = spiritGrammar.M_start;
        M_assignment                          = spiritGrammar.M_assignment;
//        M_command                             = spiritGrammar.M_command;
        M_expression                          = spiritGrammar.M_expression;
        M_compare                             = spiritGrammar.M_compare;
        M_plusMinus                           = spiritGrammar.M_plusMinus;
        M_multiplyDivide                      = spiritGrammar.M_multiplyDivide;
        M_elevate                             = spiritGrammar.M_elevate;
        M_element                             = spiritGrammar.M_element;
        M_number                              = spiritGrammar.M_number;
        M_function                            = spiritGrammar.M_function;
        M_group                               = spiritGrammar.M_group;
        M_variable                            = spiritGrammar.M_variable;
    }

    return *this;
}

// ===================================================
// Set Methods
// ===================================================
template < typename IteratorType, typename ResultsType >
inline void
ParserSpiritGrammar< IteratorType, ResultsType >::setDefaultVariables()
{
    M_variable.add( "pi" , M_PI );
    M_variable.add( "e", M_E );
}

template < typename IteratorType, typename ResultsType >
inline void
ParserSpiritGrammar< IteratorType, ResultsType >::setVariable( const std::string& name, const scalar_type& value )
{
    scalar_type *p = M_variable.find( name );
    if ( p != 0 )
        *p = value;
    else
        M_variable.add( name, value );
}

//qi::rule< IteratorType, double(double), ascii::space_type > Power, M_multiplyDivide, M_plusMinus, M_compare;
/*
        M_expression
            =   ExpressionLevel_1                          [qi::_val = qi::_1]
            >> *M_plusMinus(qi::_val)                        [qi::_val = qi::_1]
        ;

        ExpressionLevel_1
            =   M_element                                    [qi::_val = qi::_1]
            >> *M_multiplyDivide(qi::_val)                   [qi::_val = qi::_1]
        ;

        M_plusMinus
            =   qi::eps(qi::_val = qi::_r1)
            >>  (
                qi::lit('+') >> ExpressionLevel_1          [qi::_val += qi::_1]
            |   qi::lit('-') >> ExpressionLevel_1          [qi::_val -= qi::_1]
                )
        ;

        M_multiplyDivide
            =   qi::eps(qi::_val = qi::_r1)
            >>  (
                qi::lit('*') >> M_element                    [qi::_val *= qi::_1]
            |   qi::lit('/') >> M_element                    [qi::_val /= qi::_1]
                )
        ;
*/

/*#else

//! ParserSpiritGrammar - An empty implementation for boost version < 1.41
template < typename IteratorType = std::string::const_iterator, typename ResultsType = std::vector < scalar_type > >
class ParserSpiritGrammar
{
public:

    typedef IteratorType                                        iterator_Type;
    typedef boost::iterator_range< iterator_Type >              iteratorRange_Type;
    typedef ResultsType                                         results_Type;

    ParserSpiritGrammar() : M_scalar_type(0.) {}
    ParserSpiritGrammar( const ParserSpiritGrammar& ) : M_scalar_type(0.) {}
    ~ParserSpiritGrammar() {}

    ParserSpiritGrammar& operator=( const ParserSpiritGrammar& ) { return *this; }

    void clearVariables() {}

    void setDefaultVariables() {}
    void setVariable( const std::string&, const scalar_type& ) {}

    scalar_type& variable( const std::string& ) { return M_scalar_type; }

    scalar_type M_scalar_type;
};

#endif /* HAVE_BOOST_SPIRIT_QI */

} // Namespace LifeV

#endif /* Parser_SpiritGrammar_H */
