exe = alfy
nw = $(shell which noweb)
version = $(shell bash ../scripts/getVersion.sh)
date = $(shell bash ../scripts/getDate.sh)
vs = -X github.com/evolbioinf/$(exe)/util.version=$(version)
ds = -X github.com/evolbioinf/$(exe)/util.date=$(date)

all: $(exe)

$(exe): $(exe).org
	if [ "$(nw)" != "" ]; then\
		bash ../scripts/org2nw $(exe).org | notangle -R$(exe).go | gofmt > $(exe).go;\
	fi
	go build -ldflags "$(vs) $(ds)" $(exe).go
