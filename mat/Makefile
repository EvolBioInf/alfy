nw = $(shell which noweb)

all : mat

mat: mat.go
	go build mat.go
mat.go: mat.org
	if [ "$(nw)" != "" ]; then\
		awk -f ../scripts/preTangle.awk mat.org | bash ../scripts/org2nw | notangle -Rmat.go | gofmt > mat.go;\
	fi

clean:
	rm -f *.go
