/^lib/ {
    print $0, "-rpath /usr/local/lib"
}
!/^lib/
