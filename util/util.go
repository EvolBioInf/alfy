// Package util collects utility functions for alfy.
package util

import (
	"fmt"
	"github.com/evolbioinf/clio"
	"log"
	"os"
)

var version, date string

// The function Version prints the program version and other information.
func Version(name string) {
	var authors, emails, license string
	authors = "Mirjana Domazet-Loso," +
		"Beatriz Vieira Mourato," +
		"Bernhard Haubold"
	emails = "Mirjana.Domazet-Loso@fer.unizg.hr," +
		"mourato@evolbio.mpg.de," +
		"haubold@evolbio.mpg.de"
	license = "Gnu General Public License, " +
		"https://www.gnu.org/licenses/gpl.html"
	clio.PrintInfo(name, version, date,
		authors, emails, license)
	os.Exit(0)
}

// The function Check checks for errors. If yes, it prints them and exits.
func Check(e error) {
	if e != nil {
		log.Fatal(e)
	}
}

// PrepareErrorMessages takes as argument the program name and sets this as the prefix for the error messages from the log package.
func PrepareErrorMessages(name string) {
	m := fmt.Sprintf("%s - ", name)
	log.SetPrefix(m)
	log.SetFlags(0)
}
