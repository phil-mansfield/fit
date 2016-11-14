package fit

type Terminator interface {
	Stop(chains [][][][]float64) bool
}

type stepsTerminator struct {
	n int
}

func (term *stepsTerminator) Stop(chains [][][][]float64) bool {
	return len(chains[0][0]) > term.n
}

func Steps(n int) Terminator {
	return &stepsTerminator{n}
}