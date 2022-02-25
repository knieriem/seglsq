// Package seglsq is the translation of a C++ implementation
// to solve the Segmented Least Squares Problem, as described in:
//
// https://kartikkukreja.wordpress.com/2013/10/21/segmented-least-squares-problem/
//
// Original implementation:
// https://github.com/kartikkukreja/blog-codes/blob/master/src/Segmented%20Least%20Squares%20Problem.cpp
package seglsq

import (
	"math"
	"sort"
)

type Problem struct {
	points   []Point
	sorted   []Point
	segments [][]Line
}

func NewProblem(points []Point) *Problem {
	p := new(Problem)
	p.points = points
	n := len(points)

	// sort the points in non-decreasing order of x coordinate
	pts := make([]Point, n)
	copy(pts, points)
	sort.Sort(sortablePoints(pts))
	p.sorted = pts

	s := make([][]Line, n+1)
	for i := range s {
		s[i] = make([]Line, n+1)
	}
	p.segments = s
	return p
}

type Point struct {
	X int32
	Y int32
}

type sortablePoints []Point

func (p sortablePoints) Len() int           { return len(p) }
func (p sortablePoints) Less(i, j int) bool { return p[i].X < p[j].X }
func (p sortablePoints) Swap(i, j int)      { p[i], p[j] = p[j], p[i] }

type Result struct {
	Lines []Line
	Cost  float64

	x0  int32
	xn1 int32
}

// Solve tries to find an optimal solution for the problem
// with the specified cost of creating a new segment
func (p *Problem) Solve(cost float64) *Result {
	n := len(p.points)
	r := new(Result)
	r.x0 = p.points[0].X
	r.xn1 = p.points[n-1].X
	points := p.sorted
	segments := p.segments

	//
	// precompute the error terms
	//

	// cumulative[i] contains the sum of its terms for 0 <= j <= i
	cumulative := make([]terms, n+1)

	var sum terms
	for j, pt := range points {
		cNext := &cumulative[j+1]
		cNext.add(&cumulative[j], pt)

		for i := 0; i <= j; i++ {
			sum.sub(cNext, &cumulative[i])

			// seg is the segment that is best fit to points[i .. j]
			seg := &segments[i][j]
			seg.Slope, seg.Intercept = sum.lineParams(j - i + 1)

			for k := i; k <= j; k++ {
				ptk := &points[k]
				tmp := float64(ptk.Y) - seg.Slope*float64(ptk.X) - seg.Intercept
				seg.squareErr += tmp * tmp
			}
		}
	}

	//
	// Find the cost of the optimal solution
	//
	// opt[i] is the optimal solution (minimum cost) for points[0 .. i-1]
	opt := make([]float64, n+1)

	// optSeg[i] is the last segment in the optimal solution for points[0 .. i-1]
	optSeg := make([]int, n)

	for j := range points {
		k := 0
		mn := math.Inf(1)
		for i := 0; i <= j; i++ {
			tmp := segments[i][j].squareErr + opt[i]
			if tmp < mn {
				mn = tmp
				k = i
			}
		}
		opt[j+1] = mn + cost
		optSeg[j] = k
	}
	r.Cost = opt[n]

	// Find the optimal solution
	seg := make([]int, 0, n)
	for i := n - 1; i >= 0; {
		j := optSeg[i]
		seg = append(seg, i, j)
		i = j - 1
	}

	// Setup the result with optimal solution
	r.Lines = make([]Line, 0, len(seg)/2)
	for len(seg) != 0 {
		n := len(seg)
		n--
		i := seg[n]
		n--
		j := seg[n]
		seg = seg[:n]

		r.Lines = append(r.Lines, segments[i][j])
	}
	return r
}

// Intersections calculates, for the given Result,
// the intersecting points
func (r *Result) Intersections() []Point {
	isect := make([]Point, len(r.Lines)+1)

	isect[0] = r.Lines[0].point(r.x0)
	prev := r.Lines[0]
	for i, line := range r.Lines[1:] {

		isect[i+1], _ = line.intersect(&prev)
		prev = line
	}
	isect[len(isect)-1] = r.Lines[len(r.Lines)-1].point(r.xn1)
	return isect
}

// Terms are used for computing E[i][j] which is the square error of a segment
// that is best fit to points[i .. j]
type terms struct {
	x    int64
	y    int64
	xy   int64
	xSqr int64
}

// add applies a Point's coordinates to the given terms
func (sum *terms) add(t *terms, pt Point) {
	sum.x += t.x + int64(pt.X)
	sum.y += t.y + int64(pt.Y)
	sum.xy += t.xy + int64(pt.X)*int64(pt.Y)
	sum.xSqr += t.xSqr + int64(pt.X)*int64(pt.X)
}

// sub calculates the difference of min and sub terms
func (d *terms) sub(min, sub *terms) {
	d.x = min.x - sub.x
	d.y = min.y - sub.y
	d.xy = min.xy - sub.xy
	d.xSqr = min.xSqr - sub.xSqr
}

func (t *terms) lineParams(interval int) (slope, intercept float64) {
	intv := int64(interval)
	if num := intv*t.xy - t.x*t.y; num != 0 {
		denom := intv*t.xSqr - t.x*t.x
		if denom == 0 {
			slope = math.Inf(1)
		} else {
			slope = float64(num) / float64(denom)
		}
	}
	intercept = (float64(t.y) - slope*float64(t.x)) / float64(intv)
	return slope, intercept
}
