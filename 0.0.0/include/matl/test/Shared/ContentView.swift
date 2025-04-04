//
//  ContentView.swift
//  Shared
//
//  Created by Daniel Fust on 3/13/22.
//

import SwiftUI

struct ContentView: View {
    @Binding var document: testDocument

    var body: some View {
        TextEditor(text: $document.text)
    }
}

struct ContentView_Previews: PreviewProvider {
    static var previews: some View {
        ContentView(document: .constant(testDocument()))
    }
}
