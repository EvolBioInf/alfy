exe = prepAlfy
nw = $(shell which noweb)
version = $(shell bash ../scripts/getVersion.sh)
date = $(shell bash ../scripts/getDate.sh)
vs = -X github.com/evolbioinf/$(exe)/util.version=$(version)
ds = -X github.com/evolbioinf/$(exe)/util.date=$(date)

all: $(exe)

$(exe).go: $(exe).org
	if [ "$(nw)" != "" ]; then\
		bash ../scripts/org2nw $(exe).org | notangle -R$(exe).go | gofmt > $(exe).go;\
	fi
	go build -ldflags "$(vs) $(ds)" $(exe).go

$(exe)_test.go: $(exe).org
	if [ "$(nw)" != "" ]; then\
		bash ../scripts/org2nw $(exe).org | notangle -R$(exe)_test.go | gofmt > $(exe)_test.go;\
	fi

test: $(exe).go $(exe)_test.go
	@go test

clean:
	rm -f $(exe) *.go
