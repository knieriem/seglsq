package seglsq

import (
	"math"
)

type Line struct {
	Slope     float64
	Intercept float64
	squareErr float64
}

func (l *Line) intersect(l2 *Line) (i Point, ok bool) {
	defer func() {
		if x := recover(); x != nil {
			return
		}
	}()
	num := l2.Intercept - l.Intercept
	denom := l.Slope - l2.Slope
	x := num / denom
	if math.IsInf(x, 0) {
		return Point{}, false
	}
	i.X = int32(math.Round(x))
	i.Y = int32(math.Round(x*l.Slope + l.Intercept))
	return i, true
}

func (l *Line) point(x int32) Point {
	return Point{X: x, Y: l.y(x)}
}

func (l *Line) y(x int32) int32 {
	return int32(math.Round(float64(x)*l.Slope + l.Intercept))
}
