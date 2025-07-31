// Package util collects utility functions for alfy.
package util

import (
	"github.com/evolbioinf/clio"
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
